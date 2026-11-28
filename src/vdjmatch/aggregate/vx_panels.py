"""Sparse, cutoff-bound vX membership for fixed observed and mock TRB panels."""

import numpy as np
import polars as pl
from scipy.sparse import coo_matrix
from seqtree import Index, amino_acids, gapblock

from .. import db
from ..evalue.competitive import MAX_MIDPOINT_GAP_PER_CHAIN
from ..evalue.equal_tail import _eligible, _prepare
from ..evalue.vx_fixed_target import fixed_target_scores
from ..db.schema import mhc_compatible
from ..match.geometry import _gap_prior, distance_matrix, precompute_v_loop_distances
from ..match.scope import search_params

GEOMETRY_KEYS = ["cdr3b", "nativeb"]
PANEL_KEYS = ["group_id", "replicate"]


def _normalized(frame, context):
    renamed = frame.rename({"junction_aa": "cdr3b", "v": "vb"})
    return _prepare(renamed, "b", context["model"]).filter(_eligible("b", "native", 5))


class _PanelScorer:
    """Shortlist by a lower bound, then apply the frozen score to full panels."""

    def __init__(self, junctions, native_v, groups, source_log2, competitors, context, geometry,
                 cutoff, threads):
        self.junctions = list(junctions)
        self.v = list(native_v)
        self.groups = groups
        self.source_log2 = source_log2
        self.bg_seq, self.bg_v = map(list, competitors)
        self.context = context
        self.alleles, self.correction = geometry
        _, native_v_cost = context["tables"]["TRB"]
        self.weighted_v = native_v_cost + self.correction
        assert self.weighted_v.min() >= 0
        from ..cli.paired_vx_base import RULES

        self.rules = RULES["TRB"]
        self.cutoff = cutoff
        self.threads = threads
        self.ref_v = np.searchsorted(self.alleles, self.v)
        source_bound = max(weight.max() for weight in source_log2 if len(weight))
        assert self.rules["source_bonus"] >= 0
        self.radius = -cutoff + MAX_MIDPOINT_GAP_PER_CHAIN * max(0, self.rules["second_hit_gap_weight"]) + self.rules["source_bonus"] * source_bound
        self.bg_cores = [seq[3:-2] for seq in self.bg_seq]
        self.bg_group_ids = np.searchsorted(self.alleles, self.bg_v).tolist()
        self.bg_group_distances = (2 * native_v_cost).tolist()
        gap = 6 * context["distance_units"][0]
        self.bg_geometry = dict(matrix=context["count_matrix"], gap_open=gap, gap_extend=gap,
                                gap_prior=_gap_prior, threads=threads)
        # Every insertion, deletion or substitution costs at least three units
        # in the native weighted CDR3 metric. This is checked against the
        # prepared matrix, and the native single-gap objective dominates
        # unrestricted unit Levenshtein at the same minimum per-edit cost.
        aa = [amino_acids().index(letter) for letter in "ACDEFGHIKLMNPQRSTVWY"]
        c = context["costs"][np.ix_(aa, aa)]
        off_diag = c[~np.eye(len(c), dtype=bool)]
        self.edit_cost = 3 * min(off_diag[off_diag > 0].min(), context["distance_units"][0])
        assert np.min(np.diag(c)) == np.max(np.diag(c)) == 0
        assert not np.any((off_diag == 0))  # the clipped repertoire alphabet has no zero-cost substitutions
        self.budget = max(0, int(np.floor(self.radius / self.edit_cost)))
        cores = [seq[3:-2] for seq in self.junctions]
        self.nonempty = np.flatnonzero([bool(s) for s in cores])
        self.index = Index.build([cores[i] for i in self.nonempty], alphabet="aa") if self.nonempty.size else None

        self.cores = [s[3:-2] for s in self.junctions]
        self.bg_nonempty = np.flatnonzero([bool(s) for s in self.bg_cores])
        self.bg_index = Index.build([self.bg_cores[i] for i in self.bg_nonempty], alphabet="aa") if self.bg_nonempty.size else None
        self.bg_v_ids = np.asarray(self.bg_group_ids, dtype=np.intp)
        self.bg_offsets = np.asarray(self.bg_group_distances)
        refs = np.concatenate(self.groups)
        panels = np.repeat(np.arange(len(self.groups)), [len(g) for g in self.groups])
        self.membership = coo_matrix(
            (np.ones(refs.size, dtype=bool), (refs, panels)),
            shape=(len(self.cores), len(self.groups)),
        ).tocsr()

    def _distance(self, seq, v, refs, r_v):
        workers = min(self.threads, max(1, len(seq) * len(refs) // 100_000))
        return distance_matrix(seq, refs, v, r_v, species="human", threads=workers,
                               v_loop_distances=self.context, allow_short=True)

    def _pairs(self, query, reference, nonempty, index, q_v, r_v, offsets, radius):
        """Return all pairs whose unrestricted edit and V costs meet a certified bound."""
        if not query or not reference or radius < 0:
            return np.empty((0, 2), dtype=np.intp)
        q_len = np.fromiter(map(len, query), dtype=np.intp, count=len(query))
        r_len = np.fromiter(map(len, reference), dtype=np.intp, count=len(reference))
        q_nonempty = np.flatnonzero(q_len)
        gap = 6 * self.context["distance_units"][0]
        radius = np.nextafter(radius, np.inf)
        maximum = int(np.floor(2 * radius))
        edit_limit = int(np.floor(radius / self.edit_cost))
        batches = []

        def retain(q, r, bound):
            keep = bound + offsets[q_v[q], r_v[r]] <= radius
            if keep.any():
                batches.append(np.column_stack((q[keep], r[keep])))

        if index is not None and q_nonempty.size:
            # The trie uses penalty zero to mean an unrestricted weighted
            # search, so explicitly use exact sequence equality at radius zero.
            params = search_params(
                f"{edit_limit},{edit_limit},{edit_limit},{edit_limit}",
                **(dict(matrix=self.context["count_matrix"], max_penalty=maximum,
                        gap_open=gap, gap_extend=gap) if maximum else {}),
            )

            def consume(start, stop, hits):
                block = np.asarray(hits)
                if block.size:
                    retain(q_nonempty[block[:, 0]], nonempty[block[:, 1]], block[:, 2] * 0.5)

            index.search_batch_consume([query[i] for i in q_nonempty], params, consume,
                                       threads=self.threads, chunk_size=1024, max_hit_bytes=64 * 1024**2)

        # An empty core costs exactly one linear gap against a nonempty
        # core. Handle it explicitly because the trie omits empty strings.
        def empty_pairs(q, r):
            if q.size and r.size:
                qi = np.repeat(q, r.size)
                ri = np.tile(r, q.size)
                retain(qi, ri, gap * 0.5 * (q_len[qi] + r_len[ri]))

        empty_pairs(np.flatnonzero(q_len == 0), np.flatnonzero(r_len * gap <= maximum))
        empty_pairs(np.flatnonzero((q_len > 0) & (q_len * gap <= maximum)), np.flatnonzero(r_len == 0))
        return np.concatenate(batches) if batches else np.empty((0, 2), dtype=np.intp)

    def _sparse_competitor(self, query, q_v):
        radius = self.radius + self.rules["margin_per_chain"]
        pairs = self._pairs(query, self.bg_cores, self.bg_nonempty, self.bg_index, q_v,
                            self.bg_v_ids, self.bg_offsets / 2, radius)
        nearest = np.full(len(query), np.inf)
        if not pairs.size:
            return nearest
        q, r = pairs.T
        queries, group = np.unique(q, return_inverse=True)
        limit = int(np.floor(2 * np.nextafter(radius, np.inf)))
        _, counts, result = gapblock.group_audit_batch(
            [query[i] for i in queries], [self.bg_cores[i] for i in r], [limit] * len(queries),
            **self.bg_geometry,
            query_group_ids=q_v[queries].tolist(), reference_group_ids=self.bg_v_ids[r].tolist(),
            group_distances=self.bg_group_distances,
            query_partition_ids=list(range(len(queries))), reference_partition_ids=group.tolist(),
        )
        nearest[queries] = np.where(np.asarray(counts) > 0, np.asarray(result) / 2, np.inf)
        return nearest

    def batches(self, sequences, alleles):
        """Yield candidate observation indices and their panel membership masks."""
        q_v = np.searchsorted(self.alleles, alleles)
        for start in range(0, len(sequences), 32_768):
            stop = min(start + 32_768, len(sequences))
            trimmed = [s[3:-2] for s in sequences[start:stop]]
            v = q_v[start:stop]
            pairs = self._pairs(trimmed, self.cores, self.nonempty, self.index, v, self.ref_v,
                                self.weighted_v, self.radius)
            if not pairs.size:
                continue
            q, r = pairs.T
            eligible, local = np.unique(q, return_inverse=True)
            hits = coo_matrix(
                (np.ones(q.size, dtype=bool), (local, r)), shape=(eligible.size, len(self.cores))
            ).tocsr()
            panels = (hits @ self.membership).tocsc()
            background = self._sparse_competitor([trimmed[i] for i in eligible], v[eligible])
            passed = np.zeros((eligible.size, len(self.groups)), dtype=bool)
            for i, cols in enumerate(self.groups):
                rows = panels.indices[panels.indptr[i]:panels.indptr[i + 1]]
                limit = max(1, 4_000_000 // max(1, len(cols)))
                refs = [self.junctions[j] for j in cols]
                r_v = [self.v[j] for j in cols]
                for offset in range(0, rows.size, limit):
                    block = rows[offset:offset+limit]
                    positions = start + eligible[block]
                    native = self._distance([sequences[j] for j in positions], [alleles[j] for j in positions], refs, r_v)
                    delta = self.correction[np.ix_(q_v[positions], self.ref_v[cols])]
                    args = (native, delta, np.zeros(native.shape, dtype=bool), background[block], self.source_log2[i])
                    _, result = fixed_target_scores(*args, self.rules, 1)
                    passed[block, i] = result >= self.cutoff
            yield start + eligible, passed


class VxMockPanels:
    """Keep the existing observations and panel membership, and swap only the match rule."""

    def __init__(self, historical_reference, competitor_reference, epitope, mhc_a, mhc_b, mhc_class, cutoff):
        self.context = precompute_v_loop_distances(species="human")
        self.cutoff = cutoff
        original = db.load(historical_reference, species="HomoSapiens", gene="TRB", epitope=epitope)
        self.accessions = _normalized(
            original.select(pl.col("cdr3").alias("junction_aa"), "v", "reference_id"),
            self.context,
        ).group_by(GEOMETRY_KEYS).agg(
            pl.col("reference_id").drop_nulls().n_unique().alias("accessions")
        )
        raw = db.load(competitor_reference, species="HomoSapiens", gene="TRB")
        background = raw.filter(
            pl.col("epitope") != epitope, pl.col("mhc_class") == mhc_class,
            mhc_compatible("mhc_a", mhc_a), mhc_compatible("mhc_b", mhc_b),
        )
        self.background = _normalized(
            background.select(pl.col("cdr3").alias("junction_aa"), "v"), self.context
        ).select(GEOMETRY_KEYS).unique().sort(GEOMETRY_KEYS)
        if not self.background.height:
            raise ValueError("empty frozen compatible competitor pool")
        self.ready_panels = None
        self.sample_audit = []

    def initialize(self, panels, threads):
        self.inventory = panels.select(PANEL_KEYS).unique().sort(PANEL_KEYS)
        prepared = _normalized(panels, self.context).select(*PANEL_KEYS, *GEOMETRY_KEYS).unique()
        keys = prepared.select(GEOMETRY_KEYS).unique().sort(GEOMETRY_KEYS).with_row_index("key")
        prepared = prepared.join(keys, on=GEOMETRY_KEYS, validate="m:1").join(
            self.accessions, on=GEOMETRY_KEYS, how="left", validate="m:1"
        )
        groups = {
            identity: part.sort("key")
            for identity, part in prepared.partition_by(*PANEL_KEYS, as_dict=True).items()
        }
        members = []
        weights = []
        for identity in self.inventory.iter_rows():
            rows = groups.get(identity, prepared.head(0))
            bases = rows["accessions"].fill_null(1).clip(lower_bound=1).to_numpy()
            members.append(rows["key"].to_numpy())
            weights.append(np.log2(bases) if identity[1] == -1 else np.zeros(len(rows)))
        from ..cli.paired_vx_base import _regional_corrections

        if not keys.height:
            raise ValueError("no resolved receptor key in any vX reference panel")
        geometry = _regional_corrections(self.context, "TRB")
        self.scorer = _PanelScorer(keys["cdr3b"], keys["nativeb"], members, weights,
                                  (self.background["cdr3b"], self.background["nativeb"]),
                                  self.context, geometry, self.cutoff, threads)
        self.ready_panels = panels
        print(f"{len(keys):,} unique panel keys; {len(self.background):,} compatible competitor keys; "
              f"edit shortlist bound {self.scorer.budget}", flush=True)

    def count(self, observations, panels, *, threads):
        if self.ready_panels is None:
            self.initialize(panels, threads)
        if panels is not self.ready_panels:
            raise ValueError("vX requires the same fixed panels throughout the run")
        queries = _normalized(observations, self.context).group_by(GEOMETRY_KEYS).agg(
            pl.len().alias("source_rows"),
            pl.col("count").sum().alias("reads"),
            pl.col("count").alias("abundances"),
        ).sort(GEOMETRY_KEYS)
        query_seq = queries["cdr3b"].to_list()
        query_v = queries["nativeb"].to_list()
        abundances = [[] for _ in range(self.inventory.height)]
        shortlisted = 0
        for positions, matches in self.scorer.batches(query_seq, query_v):
            shortlisted += len(positions)
            for j in np.flatnonzero(matches.any(axis=0)):
                selected = positions[matches[:, j]]
                abundance = queries[selected.tolist()]["abundances"].explode().to_list()
                abundances[j].extend(abundance)
        result = self.inventory.with_columns(
            pl.Series("matched_records", [len(v) for v in abundances], dtype=pl.Int64),
            pl.Series("matched_reads", [sum(v) for v in abundances], dtype=pl.Int64),
            pl.Series("matched_abundances", abundances, dtype=pl.List(pl.Int64)),
        )
        self.sample_audit.append({
            "source_rows": observations.height,
            "resolved_rows": int(queries["source_rows"].sum()),
            "distinct_geometry_keys": queries.height,
            "shortlisted_geometry_keys": shortlisted,
        })
        print(f"Geometry keys {queries.height:,}; shortlisted {shortlisted:,} across all panels", flush=True)
        return result
