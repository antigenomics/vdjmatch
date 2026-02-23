#!/usr/bin/env python
from __future__ import annotations

import argparse
import os
import tempfile
from typing import List, Optional, Tuple

import numpy as np
import pandas as pd
from tqdm import tqdm

from tcrtrie import Trie
from mir.basic.pgen import OlgaModel


AA20 = set("ACDEFGHIKLMNPQRSTVWY")


def strip_allele(gene_call: Optional[str]) -> str:
    if not gene_call:
        return ""
    s = str(gene_call).strip()
    if "*" in s:
        s = s.split("*", 1)[0]
    return s.strip()


def is_valid_cdr3(s: str) -> bool:
    if not s:
        return False
    s = s.strip()
    if not (6 <= len(s) <= 30):
        return False
    return all(ch in AA20 for ch in s)


def stable_group_ids(df: pd.DataFrame) -> Tuple[pd.Series, List[Tuple[str, str, str]]]:
    uniq = df[["species", "chain", "group"]].drop_duplicates().reset_index(drop=True)

    gid_map = {}
    gid_to_key = []

    for i, row in uniq.iterrows():
        key = (row["species"], row["chain"], row["group"])
        gid_map[key] = i
        gid_to_key.append(key)

    gids = df.apply(
        lambda r: gid_map[(r["species"], r["chain"], r["group"])],
        axis=1,
    )

    return gids, gid_to_key


def main():
    ap = argparse.ArgumentParser()

    ap.add_argument("input_tsv")
    ap.add_argument("-n", type=int, required=True)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--batch-size", type=int, default=10000)
    ap.add_argument("--gene", default="TRB")
    ap.add_argument("--group-col", default="antigen.epitope")
    ap.add_argument("--junction-col", default="junction_aa")
    ap.add_argument("--v-col", default="v_call")
    ap.add_argument("--j-col", default="j_call")
    ap.add_argument("--species-col", default="species")
    ap.add_argument("--chain-col", default="gene")
    ap.add_argument("--match-v", action="store_true")
    ap.add_argument("--match-j", action="store_true")
    ap.add_argument("--out", default="background_matches.tsv")
    ap.add_argument("--max-sub", type=int, default=1)
    ap.add_argument("--max-ins", type=int, default=0)
    ap.add_argument("--max-del", type=int, default=0)
    ap.add_argument("--max-edits", type=int, default=None)

    args = ap.parse_args()

    df = pd.read_csv(args.input_tsv, sep="\t", dtype=str, keep_default_na=False)

    df = df[
        [
            args.junction_col,
            args.v_col,
            args.j_col,
            args.group_col,
            args.species_col,
            args.chain_col,
        ]
    ].copy()

    df.columns = ["junction_aa", "v_call", "j_call", "group", "species", "chain"]

    for c in df.columns:
        df[c] = df[c].astype(str).str.strip()

    df["v_call"] = df["v_call"].map(strip_allele)
    df["j_call"] = df["j_call"].map(strip_allele)

    df = df[df["junction_aa"].map(is_valid_cdr3)]

    if df.empty:
        raise RuntimeError("Input empty after filtering")

    gids, gid_to_key = stable_group_ids(df)
    df["__group_id"] = gids.astype(int)

    with tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False) as tmp:
        tmp_path = tmp.name
        df[["junction_aa", "v_call", "j_call", "__group_id"]].to_csv(
            tmp_path, sep="\t", index=False
        )

    try:
        trie = Trie(tmp_path)

        model = OlgaModel(chain=args.gene)
        np.random.seed(args.seed)

        counts = [0] * len(gid_to_key)

        batch_q = []
        batch_v = []
        batch_j = []

        pbar = tqdm(total=args.n, unit="seq")

        produced = 0

        while produced < args.n:
            k = min(args.batch_size, args.n - produced)
            samples = model.generate_sequences_with_meta(k, False)

            for rec in samples:
                cdr3 = str(rec.get("cdr3", "")).strip()
                if not is_valid_cdr3(cdr3):
                    continue

                v = strip_allele(rec.get("v_gene", ""))
                j = strip_allele(rec.get("j_gene", ""))

                batch_q.append(cdr3)

                if args.match_v:
                    batch_v.append(v)
                if args.match_j:
                    batch_j.append(j)

                produced += 1
                if produced >= args.n:
                    break

            if len(batch_q) >= args.batch_size or produced >= args.n:

                v_filters = batch_v if args.match_v else None
                j_filters = batch_j if args.match_j else None

                results = trie.SearchGroupIdsForAll(
                    batch_q,
                    maxSubstitution=args.max_sub,
                    maxInsertion=args.max_ins,
                    maxDeletion=args.max_del,
                    maxEdits=args.max_edits,
                    vGeneFilters=v_filters,
                    jGeneFilters=j_filters,
                    unique=True,
                )

                for gids_list in results:
                    for gid in gids_list:
                        counts[int(gid)] += 1

                pbar.update(len(batch_q))

                batch_q.clear()
                batch_v.clear()
                batch_j.clear()

        pbar.close()

        out_rows = []
        for gid, (sp, ch, grp) in enumerate(gid_to_key):
            m = counts[gid]
            out_rows.append(
                {
                    "species": sp,
                    "chain": ch,
                    args.group_col: grp,
                    "n_background": args.n,
                    "matched_background": m,
                    "match_fraction": m / args.n,
                }
            )

        pd.DataFrame(out_rows).to_csv(args.out, sep="\t", index=False)

    finally:
        try:
            os.remove(tmp_path)
        except OSError:
            pass


if __name__ == "__main__":
    main()