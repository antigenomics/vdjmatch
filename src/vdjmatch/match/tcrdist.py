"""TCRdist3 default human/mouse alpha/beta distance using native gap-block scoring.

Reproduces 3*CDR3 + aligned CDR1/CDR2/CDR2.5, including TCRdist3's clipped
BLOSUM matrix and restricted gap positions. This distance is not calibrated evidence.
"""
from __future__ import annotations

import csv
from importlib import resources

import numpy as np
from seqtree import SubstitutionMatrix, amino_acids, gapblock


def kneedle_radius(distances):
    """Strongest confirmed knee of sorted nearest distances, or None.

    Offline Kneedle on an increasing distance curve: unit-square normalization,
    convex-curve difference x-y, local maxima and fixed sensitivity S=1. The
    observed monotone curve is used without smoothing or fitted parameters.
    Multiplying all distances by c>0 multiplies the returned radius by c.
    Flat/linear/short curves have no boundary; endpoints cannot be knees.
    """
    from scipy.signal import find_peaks

    values = np.asarray(distances, dtype=float)
    if values.ndim != 1 or not np.isfinite(values).all() or np.any(values < 0):
        raise ValueError('nearest distances must be a finite nonnegative vector')
    if len(values) < 3 or np.ptp(values) == 0:
        return None
    values = np.sort(values)
    difference = np.linspace(0, 1, len(values)) - (values-values[0])/np.ptp(values)
    peaks = find_peaks(difference)[0]
    minima = find_peaks(-difference)[0]
    step = 1/(len(values)-1)
    for peak in sorted(peaks, key=lambda i: difference[i], reverse=True):
        if difference[peak] <= 8*np.finfo(float).eps:
            continue
        following = minima[minima > peak]
        end = int(following[0]) if len(following) else len(values)-1
        if np.any(difference[peak+1:end+1] < difference[peak]-step):
            return float(values[peak])
    return None


def regularize_reference_radii(radii, nearest_distances, *, factor=.5):
    """Floor knees within one species/locus/geometry; leave absent knees absent.

    Use factor times the median positive knee, falling back to positive nearest
    distances when all knees are zero. Both scales transform with distance units.
    Without any positive scale a zero knee is unavailable, never an exact-only ball.
    """
    raw = np.asarray(radii, dtype=float)
    distances = np.asarray(nearest_distances, dtype=float)
    if (raw.ndim != 1 or distances.ndim != 1 or np.isinf(raw).any()
            or np.isinf(distances).any() or np.any(raw < 0) or np.any(distances < 0)
            or not np.isfinite(factor) or factor <= 0):
        raise ValueError('nonnegative vectors and a finite positive knee floor factor required')
    positive = raw[np.isfinite(raw) & (raw > 0)]
    source = 'positive_knee_median'
    if not len(positive):
        positive = distances[np.isfinite(distances) & (distances > 0)]
        source = 'positive_nearest_median'
    if not len(positive):
        return np.full_like(raw, np.nan), None, 'no_positive_scale'
    floor = float(factor*np.median(positive))
    if not np.isfinite(floor) or floor <= 0:
        raise ValueError('knee floor must be finite and positive')
    return np.maximum(raw, floor), floor, source


def reference_structure(distances, *, radius_factor=22.5, costs=None):
    """Reference neighbour coverage and normalized convex bend, without labels.

    Coverage counts all reference keys, including keys without a nonexact
    neighbour. The radius is a factor of the shared substitution-matrix scale.
    The bend is max(x-y) of the sorted finite distances normalized to a unit
    square; flat/short curves have zero bend. Both features lie in [0,1].
    """
    values = np.asarray(distances, dtype=float)
    if values.ndim != 1 or np.isinf(values).any() or np.any(values < 0):
        raise ValueError('nonnegative nearest distances or NaN required')
    if not np.isfinite(radius_factor) or radius_factor <= 0:
        raise ValueError('positive finite matrix-relative radius factor required')
    radius = radius_factor*distance_units(costs)[0]
    coverage = float(np.count_nonzero(values <= radius)/len(values)) if len(values) else 0.
    finite = np.sort(values[np.isfinite(values)])
    shape = 0.
    if len(finite) >= 3 and np.ptp(finite) > 0:
        shape = float(np.max(np.linspace(0, 1, len(finite))-(finite-finite[0])/np.ptp(finite)))
    return coverage, shape


def reference_nearest_distances(frame, *, paired=False, species='human', threads=1, ctrim=3, alphabet='full', stem_depth=0, stem_weight=1.0, matrix='tcrdist',v_loop_distances=None, v_regions=None):
    """Nonexact same-population nearest distances from one native batch.

    Input columns are seq0/v0, plus seq1/v1 for linked TRA/TRB. Full junctions
    define self-exclusion; a pair excludes only the same full linked receptor.
    NaN means no nonexact reference. No K cap or dense receptor distance matrix
    is used. The paired reducer uses the existing sum against one reference pair.
    """
    v_loop_distances = _v_region_context(v_regions, v_loop_distances, species, alphabet, matrix)
    import polars as pl

    costs = _substitution_costs(alphabet,matrix)
    if not frame.height:
        return np.empty(0, dtype=float)
    axes = range(2 if paired else 1)
    if v_loop_distances is not None and v_loop_distances.get('species') != species:
        raise ValueError('precomputed V loops require matching species')
    model = load_v_loops(species) if v_loop_distances is None else v_loop_distances['model']
    scale,quantum=distance_units(costs)
    bound = 0
    for i in axes:
        genes = resolve_v_alleles(frame[f'v{i}'], model,
                                 locus=('TRA', 'TRB')[i] if paired else None)
        if genes.null_count():
            raise ValueError('reference V calls must resolve to the selected locus')
        width = max(len(model[g]['cdr1']+model[g]['cdr2']+model[g]['cdr25'])
                    for g in genes.unique())
        bound += int(np.ceil(stem_weight*scale)) if stem_depth else 0
        bound += scale*(3*frame[f'seq{i}'].str.len_chars().max()+width)
    sequence_columns = [f'seq{i}' for i in axes]
    sizes = frame.group_by(sequence_columns).len().rename({'len': '_self_count'})
    expected = frame.join(sizes, on=sequence_columns, how='left', validate='m:1',
                          maintain_order='left').select(frame.height-pl.col('_self_count')).to_series().to_numpy()
    if not paired:
        q, r, options = _total_distance_options(frame['seq0'], frame['seq0'],
            frame['v0'], frame['v0'], species, threads, ctrim, 'uniform', alphabet=alphabet, stem_depth=stem_depth, stem_weight=stem_weight, matrix=matrix,v_loop_distances=v_loop_distances)
        _, counts, nearest = gapblock.group_audit_batch(q, r, [200*bound]*len(q),
                                                       exclude_exact=True, **options)
        if not np.array_equal(counts, expected):
            raise ArithmeticError('reference nearest bound did not cover the population')
        result = np.asarray(nearest, dtype=float)/200
        result[expected == 0] = np.nan
        return result
    cutoffs=list(range(0,200*bound+1,200*quantum))
    if frame.height*len(cutoffs)*8 > 256*1024**2:
        raise ValueError('paired reference count curve exceeds256MiB; select a smaller target')
    counts, _ = paired_total_distance_mass_batch(
        frame['seq0'], frame['seq1'], frame['seq0'], frame['seq1'],
        frame['v0'], frame['v1'], frame['v0'], frame['v1'],
        [cutoffs]*frame.height,
        species=species, threads=threads, ctrim=ctrim, exclude_exact=True, alphabet=alphabet, stem_depth=stem_depth, stem_weight=stem_weight, matrix=matrix,v_loop_distances=v_loop_distances)
    curve = np.asarray(counts, dtype=np.uint64)
    if not np.array_equal(curve[:, -1], expected):
        raise ArithmeticError('paired reference nearest bound did not cover the population')
    result = quantum*np.argmax(curve > 0, axis=1).astype(float)
    result[expected == 0] = np.nan
    return result


def _gap_prior(i, gap, width):
    return _gap_window(i, gap, width, 5)


def _gap_prior_trim3(i, gap, width):
    return _gap_window(i, gap, width, 6)


def _gap_window(i, gap, width, trimmed, ntrim=3):
    if not gap:
        return 0
    short = width-gap+trimmed  # restore full junction coordinates
    low, high = 5, short-5
    while low > high:
        low -= 1
        high += 1
    return 0 if low-ntrim <= i <= high-ntrim else gapblock.UNREACHABLE


def _full_gap_prior(i, gap, width):
    return _gap_window(i, gap, width, 0, ntrim=0)


def load_v_loops(species='human'):
    """Bundled, source-pinned aligned loop templates keyed by exact allele."""
    if species not in ('human','mouse'):
        raise ValueError('TCRdist species must be human or mouse')
    path=resources.files('vdjmatch.resources')/'tcrdist'/f'{species}_v_loops.tsv'
    with path.open() as f:
        return {r['allele']:r for r in csv.DictReader(f,delimiter='\t')}


def resolve_v_alleles(calls, model, *, locus=None):
    """Resolve exact/model-proven co-locus spellings, preserving callers' raw calls.

    Omitted allele means *01. Only a unique modeled co-locus gene can supply an
    omitted /DV suffix or slash; no gene-family prefix or first comma call is chosen.
    Unknown alleles, missing inputs, cross-locus calls and ambiguous aliases are null.
    Mouse -DV and /DV spellings merge only when their aligned templates agree.
    """
    import polars as pl
    aliases,exact_targets=_v_allele_aliases(model)
    values=pl.Series(calls,dtype=pl.String)
    unique=values.drop_nulls().unique()
    resolved=[]
    for raw in unique:
        token=raw.strip()
        exact=token if '*' in token else token+'*01'
        candidates={exact_targets[exact]} if exact in model else aliases.get(token,set())
        target=next(iter(candidates)) if len(candidates)==1 else None
        resolved.append(target if target is not None and (locus is None or model[target]['locus']==locus) else None)
    return values.replace_strict(unique,pl.Series(resolved,dtype=pl.String),default=None,return_dtype=pl.String)


def _v_allele_aliases(model):
    """Model-proven spelling map; independent of query/reference populations."""
    aliases={};exact_targets={}
    for key,meta in model.items():
        gene,allele=key.split('*')
        canonical=gene.replace('-DV','/DV')+'*'+allele
        if canonical not in model or any(meta[c]!=model[canonical][c] for c in ('locus','cdr1','cdr2','cdr25')):
            canonical=key
        exact_targets[key]=canonical
        forms={gene}
        if '/DV' in canonical:
            full=canonical.split('*')[0]
            forms.update([full,full.replace('/',''),full.split('/DV')[0]])
        for form in forms:
            tokens=[form+'*'+allele]+([form] if allele=='01' else [])
            for token in tokens:
                aliases.setdefault(token,set()).add(canonical)
    return aliases,exact_targets


def _v_allele_lookup(model):
    aliases,exact_targets=_v_allele_aliases(model)
    tokens=set(aliases)|set(exact_targets)
    lookup={}
    for token in tokens:
        exact=token if '*' in token else token+'*01'
        candidates={exact_targets[exact]} if exact in model else aliases.get(token,set())
        lookup[token]=next(iter(candidates)) if len(candidates)==1 else None
    return lookup


def _resolved_v_calls(calls,model,v_loop_distances=None,*,locus=None):
    if v_loop_distances is None or 'resolved_calls' not in v_loop_distances:
        return resolve_v_alleles(calls,model,locus=locus).to_list()
    lookup=v_loop_distances['resolved_calls']
    genes=[lookup.get(v.strip()) if isinstance(v,str) else None for v in calls]
    if locus is not None:
        genes=[g if g is not None and model[g]['locus']==locus else None for g in genes]
    return genes


def _prepared_costs(alphabet,matrix,v_loop_distances):
    if (v_loop_distances is not None and v_loop_distances.get('alphabet')==alphabet
            and v_loop_distances.get('matrix')==matrix):
        return v_loop_distances['costs']
    return _substitution_costs(alphabet,matrix)


def _substitution_costs(alphabet='full', matrix='tcrdist'):
    """Original costs or their ordered IMGT-class compatibility projection.

    Keep raw sequence encodings: zero class cost does not imply raw identity.
    The same matrix supplies CDR3, V loops, gaps and radius/bin units.
    """
    if alphabet not in ('full', 'imgt'):
        raise ValueError('alphabet must be full or imgt')
    if matrix not in ('tcrdist', 'blosum62', 'pam250', 'mj', 'vdjam'):
        raise ValueError('matrix must be tcrdist, blosum62, pam250, mj or vdjam')
    if matrix != 'tcrdist' and alphabet != 'full':
        raise ValueError('unclipped matrix requires the full alphabet')
    symbols=amino_acids()
    if matrix == 'vdjam':
        from .scoring import load_vdjam
        b=load_vdjam()
    else:
        b=getattr(SubstitutionMatrix, {'tcrdist':'blosum62', 'mj':'structural'}.get(matrix,matrix))()
    if matrix != 'tcrdist':
        return np.asarray([[b.penalty(a,c) for c in symbols] for a in symbols],dtype=np.int32)
    costs=np.asarray([[0 if a==c else max(0,min(4,4-b.similarity(a,c)))
                       for c in symbols] for a in symbols],dtype=np.int32)
    if alphabet=='full':return costs
    from ..evalue.epitope_rules import IMGT_CLASSES
    groups={}
    for i,a in enumerate(symbols):groups.setdefault(IMGT_CLASSES.get(a,a),[]).append(i)
    result=costs.copy()
    for a,ai in groups.items():
        for b,bi in groups.items():
            block=np.ix_(ai,bi)
            result[block]=0 if a==b else costs[block].min()
    return result


def distance_units(costs=None):
    """Maximum cost and integer quantum of the shared CDR3/V-loop matrix.

    Radius and bin factors use this common distance unit. The default clipped
    BLOSUM62 unit is 4; the ordinary CDR3 gap cost is three such units.
    """
    costs=_substitution_costs() if costs is None else costs
    scale=np.max(costs)
    if not np.isfinite(costs).all() or np.any(costs<0) or np.any(costs!=np.floor(costs)) or scale<=0:
        raise ValueError('substitution matrix requires a finite positive distance scale')
    return int(scale),int(np.gcd.reduce(costs[costs>0].astype(np.int64)))


def distance_cutoff(radius):
    """Floor a decimal distance to native units, without repeated float rounding."""
    from decimal import Decimal, ROUND_FLOOR
    if not np.isfinite(radius) or radius<0:
        raise ValueError('radius must be finite and nonnegative')
    return int((Decimal(str(radius))*200).to_integral_value(rounding=ROUND_FLOOR))


def _cdr3_options(threads=1, costs=None, ctrim=2, v_loop_distances=None):
    """Native CDR3 geometry; native scores equal2/3 of weighted TCRdist CDR3."""
    costs=_substitution_costs() if costs is None else costs
    matrix=(v_loop_distances['cdr3_matrix'] if v_loop_distances is not None and costs is v_loop_distances['costs']
            else SubstitutionMatrix.from_similarity((-costs).tolist()))
    return dict(matrix=matrix,
                gap_open=2*distance_units(costs)[0],gap_extend=2*distance_units(costs)[0],
                gap_prior=_gap_prior if ctrim==2 else _gap_prior_trim3,threads=threads)


def _distance_inputs(queries, references, query_v, reference_v, species, threads, ctrim, *, dense=False, v_loop_distances=None):
    """Shared junction/gene availability contract for dense distances and counts."""
    if isinstance(ctrim,bool) or not isinstance(ctrim,int) or ctrim not in (2,3):
        raise ValueError('ctrim must be 2 (TCRdist3) or 3 (symmetric trim)')
    if isinstance(threads,bool) or not isinstance(threads,int) or threads<1:
        raise ValueError('threads must be positive')
    queries, references, query_v, reference_v = map(list,(queries,references,query_v,reference_v))
    if len(queries)!=len(query_v) or len(references)!=len(reference_v):
        raise ValueError('each junction requires a V call')
    if dense and 4*len(queries)*len(references)>64*1024**2:
        raise ValueError('distance output exceeds64MiB; use contiguous query batches')
    canonical=set('ACDEFGHIKLMNPQRSTVWY')
    if any(not isinstance(s,str) or len(s)<8 or not set(s)<=canonical for s in queries+references):
        raise ValueError('TCRdist junctions require canonical amino acids and length>=8')
    if v_loop_distances is not None and v_loop_distances.get('species') != species:
        raise ValueError('precomputed V loops require matching species')
    model=load_v_loops(species) if v_loop_distances is None else v_loop_distances['model']
    genes=_resolved_v_calls(query_v+reference_v,model,v_loop_distances)
    if any(v not in model for v in genes):
        raise ValueError(f'unknown {species} V allele; handle missing gene evidence before scoring')
    if len({model[v]['locus'] for v in genes})>1:
        raise ValueError('distance batch must contain one locus')
    return queries, references, model, genes


def precompute_v_loop_distances(*, species='human', alphabet='full', matrix='tcrdist', v_regions='cdr12plus'):
    """Exact aligned V-loop costs for reuse within one annotation invocation.

    The small allele tables are immutable inputs to later population batches;
    no sample predictions, process-global state or persisted cache are retained.
    TRA/TRB have separate tables because their aligned widths can differ.
    ``v_regions="cdr12"`` omits CDR2.5; the default ``cdr12plus`` retains it.
    """
    if v_regions not in ('cdr12plus', 'cdr12'):
        raise ValueError('v_regions must be cdr12plus or cdr12')
    model=load_v_loops(species)
    selected_model = (model if v_regions == 'cdr12plus' else
                      {v: dict(meta, cdr25='') for v, meta in model.items()})
    costs=_substitution_costs(alphabet,matrix)
    tables={}
    for locus in ('TRA','TRB'):
        alleles=np.asarray(sorted(v for v in model if model[v]['locus']==locus))
        distances,_=_v_loop_groups(alleles.tolist(),selected_model,costs)
        distances.setflags(write=False)
        alleles.setflags(write=False)
        tables[locus]=(alleles,distances)
    costs.setflags(write=False)
    return {'species':species,'alphabet':alphabet,'matrix':matrix,'model':model,'v_regions':v_regions,
            'costs':costs,'tables':tables,'resolved_calls':_v_allele_lookup(model),
            'distance_units':distance_units(costs),
            'count_matrix':SubstitutionMatrix.from_similarity((-3*costs).tolist()),
            'cdr3_matrix':SubstitutionMatrix.from_similarity((-costs).tolist())}


def _v_region_context(v_regions, context, species, alphabet, matrix):
    """Resolve optional aligned-region choice without changing default preparation.

    None inherits a supplied context, otherwise selects CDR1+CDR2+CDR2.5.
    ``cdr12`` omits only CDR2.5; allele resolution and junction geometry stay common.
    """
    if v_regions is None:
        return context
    if v_regions not in ('cdr12plus', 'cdr12'):
        raise ValueError('v_regions must be cdr12plus or cdr12')
    if context is not None:
        if context.get('v_regions', 'cdr12plus') != v_regions:
            raise ValueError('precomputed V loops require matching V regions')
        return context
    if v_regions == 'cdr12plus':
        return None
    return precompute_v_loop_distances(species=species, alphabet=alphabet,
                                      matrix=matrix, v_regions=v_regions)


def _v_loop_groups(genes, model, d, v_loop_distances=None):
    """Small aligned-loop table plus per-input group IDs, using modeled alleles."""
    if not genes:
        return np.zeros((1,1),dtype=np.int32), np.empty(0,dtype=np.int32)
    if v_loop_distances is not None:
        if not np.array_equal(d,v_loop_distances['costs']):
            raise ValueError('precomputed V loops require matching substitution costs')
        alleles,indices=np.unique(genes,return_inverse=True)
        locus=model[alleles[0]]['locus']
        all_alleles,distances=v_loop_distances['tables'][locus]
        selected=np.searchsorted(all_alleles,alleles)
        if np.any(selected>=len(all_alleles)) or not np.array_equal(all_alleles[selected],alleles):
            raise ValueError('precomputed V loops do not cover the requested alleles')
        return distances[np.ix_(selected,selected)],indices
    # TCRdist's unknown/dot/star row is zero. Score unique allele pairs once per call.
    canonical=set('ACDEFGHIKLMNPQRSTVWY')
    alphabet=amino_acids()
    lut=np.zeros((128,128),dtype=np.int32)
    for a in canonical:
        for c in canonical:lut[ord(a),ord(c)]=d[alphabet.index(a),alphabet.index(c)]
    alleles,indices=np.unique(genes,return_inverse=True)
    loops=[model[v]['cdr1']+model[v]['cdr2']+model[v]['cdr25'] for v in alleles]
    width=len(loops[0])
    codes=np.frombuffer(''.join(loops).encode('ascii'),dtype=np.uint8).reshape(len(alleles),width)
    loop_dist=np.zeros((len(alleles),len(alleles)),dtype=np.int32)
    for pos in range(width):loop_dist+=lut[codes[:,None,pos],codes[None,:,pos]]
    return loop_dist, indices


STEM_DEPTH_DEFAULT = 4


def stem_penalty(queries, references, *, depth=STEM_DEPTH_DEFAULT, weight=1.0, costs=None):
    """Joint V-J hairpin stem cost, in the shared substitution-matrix scale.

    The CDR3 loop folds as a hairpin, so junction position ``i`` sits opposite position
    ``L-1-i``. A *stem level* is that pair of positions, and a level is only conserved
    when **both** of its residues are conserved, so the two residue costs combine as a
    maximum rather than a sum: one substitution anywhere in the pair breaks the level.

        d += weight * (1/depth) * SUM over i < depth of
             max( cost(q[i], r[i]), cost(q[Lq-1-i], r[Lr-1-i]) )

    The level costs come from the same substitution table the junction alignment
    and the V loops use, so the term carries no private unit and no hand-written
    distance: ``weight=1`` means one mean level cost per pair, which spans the same
    0..``distance_units(costs)[0]`` range as one aligned residue. It is deliberately NOT
    rescaled by ``Pmax`` -- that is what the measured ``stem4`` arm does, and rescaling
    changes the recorded geometry.
    ``depth`` is a declared residue count, not a per-clonotype boundary; ``depth=0``
    disables the term and returns zeros. Junctions are used untrimmed: the stem is
    anchored on the conserved ends that the alignment trim removes, which is why it adds
    information the gap block cannot see.

    Returns an integer ``(len(queries), len(references))`` matrix in native cost units,
    ready to add to :func:`distance_matrix` output.

    Naming: ``bench/holdout_stem_apex.py`` uses "stem" for a different object -- the
    k-mers touching the first/last three residues, contrasted there against the loop
    "apex". The two are not the same construction, but they do cover the same germline
    frame, and that benchmark's finding is consistent with this term's behaviour: the
    stem region tracks the V marginal and is robust rather than sharp, which is why this
    term gives small, consistent gains and helps most on TRA, where the junction is
    80.5 per cent germline against 63.5 per cent on TRB.
    """
    if isinstance(depth, bool) or not isinstance(depth, int) or depth < 0:
        raise ValueError('stem depth must be a nonnegative integer')
    if isinstance(weight, bool) or not isinstance(weight, (int, float)) or not np.isfinite(weight) or weight < 0:
        raise ValueError('stem weight must be finite and nonnegative')
    queries, references = list(queries), list(references)
    out = np.zeros((len(queries), len(references)), dtype=np.int32)
    if not depth or not weight or not queries or not references:
        return out
    costs = _substitution_costs() if costs is None else costs
    if any(len(s) < 2 * depth for s in queries + references):
        raise ValueError(f'stem depth {depth} requires junctions of length >= {2 * depth}')
    alphabet = amino_acids()
    index = {c: i for i, c in enumerate(alphabet)}
    def ends(seqs):
        head = np.array([[index[s[i]] for i in range(depth)] for s in seqs], dtype=np.int64)
        tail = np.array([[index[s[len(s) - 1 - i]] for i in range(depth)] for s in seqs],
                        dtype=np.int64)
        return head, tail
    qh, qt = ends(queries)
    rh, rt = ends(references)
    total = np.zeros((len(queries), len(references)), dtype=np.int64)
    for level in range(depth):
        total += np.maximum(costs[qh[:, level][:, None], rh[None, :, level]],
                            costs[qt[:, level][:, None], rt[None, :, level]])
    result = np.rint(weight * total / depth)
    if np.any(result > np.iinfo(np.int32).max // 200):
        raise ValueError('stem weight exceeds native distance units')
    return result.astype(np.int32)


def distance_matrix(queries, references, query_v, reference_v, *, species='human', threads=1, return_cdr3=False, ctrim=2, v_loop_weight=1, alphabet='full', stem_depth=0, stem_weight=1.0, matrix='tcrdist', v_loop_distances=None, v_regions=None):
    """TCRdist3 default distances in a bounded native batch, maximum64MiB output.

    Exact alleles are used; a missing allele suffix explicitly means *01, as in the
    comparator. Model-proven co-locus aliases resolve without changing raw calls.
    Unknown alleles, ambiguous names, cross-locus inputs and junctions shorter than8 fail
    instead of silently losing gene evidence. Dots in aligned loops carry the original
    TCRdist3 unknown-symbol cost zero (also for stars). Human templates are28 aligned
    positions; mouse TRA templates are29. Callers handle missing calls explicitly.
    ``return_cdr3=True`` returns total and weighted CDR3 matrices separately.
    ``v_loop_weight=0`` is the fixed CDR3-only ablation; eligibility stays identical.
    ``v_regions="cdr12"`` retains aligned CDR1/CDR2 only. None inherits a
    supplied prepared context, otherwise keeps CDR1/CDR2/CDR2.5.
    ``ctrim=3`` tests removal of three residues at both junction ends; the default
    remains the TCRdist3 three/two trim. Full junction identity is handled by callers.
    ``stem_depth`` adds the joint V-J hairpin stem term of :func:`stem_penalty` on the
    untrimmed junctions; 0 (the default) reproduces TCRdist3 exactly. The returned
    ``cdr3`` matrix under ``return_cdr3`` stays the junction-only distance and never
    carries the stem, which is a joint V-J term rather than a CDR3 alignment cost.
    """
    v_loop_distances = _v_region_context(v_regions, v_loop_distances, species, alphabet, matrix)
    if type(v_loop_weight) is not int or v_loop_weight not in (0,1):
        raise ValueError("v_loop_weight must be integer0 or1")
    queries, references, model, genes = _distance_inputs(
        queries, references, query_v, reference_v, species, threads, ctrim, dense=True, v_loop_distances=v_loop_distances)
    d=_prepared_costs(alphabet,matrix,v_loop_distances)
    stem_penalty([],[],depth=stem_depth,weight=stem_weight,costs=d)
    if not queries or not references:
        out=np.empty((len(queries),len(references)),dtype=np.int32)
        return (out,out.copy()) if return_cdr3 else out
    # Integer similarities -d produce exactly2*d under seqtree's Gram transform.
    cdr3=np.array(gapblock.score_matrix([s[3:-ctrim] for s in queries],[s[3:-ctrim] for s in references],
        **_cdr3_options(threads,d,ctrim,v_loop_distances)),copy=True)
    cdr3//=2
    cdr3*=3
    out=cdr3.copy() if return_cdr3 else cdr3
    loop_dist, indices = _v_loop_groups(genes, model, d, v_loop_distances)
    out+=v_loop_weight*loop_dist[indices[:len(queries),None],indices[None,len(queries):]]
    if stem_depth:
        out+=stem_penalty(queries,references,depth=stem_depth,weight=stem_weight,costs=d)
    return (out,cdr3) if return_cdr3 else out


def total_distance_count_batch(queries, references, query_v, reference_v, thresholds, *,
                               species='human', threads=1, ctrim=2,
                               position_weighting='uniform', exclude_exact=False, linear_mass=False,
                               query_boundaries=None, reference_boundaries=None, v_loop_weight=1,
                               same_v_gene=False, alphabet='full', stem_depth=0, stem_weight=1.0, matrix='tcrdist', shared_thresholds=False, v_loop_distances=None, v_regions=None):
    """Count full-junction/V-loop balls without a dense query/reference matrix.

    ``thresholds`` contains one integer list per query, in native units200 times
    the ordinary TCRdist distance. Results preserve threshold order/duplicates.
    ``shared_thresholds=True`` accepts one flat cutoff list for the whole batch.
    ``uniform`` reproduces distance_matrix with three N-terminal and ctrim (2/3)
    C-terminal positions zero-weighted. ``significance`` instead uses the shipped
    full end-anchored profile, quantized max(1,round(100*w)), including flanks;
    it rejects ctrim3 rather than combining trimming and decay. Both retain the
    full-coordinate TCRdist moving gap window and linear gap cost12 per residue.

    The clipped substitution matrix contributes threefold CDR3 distance; aligned
    CDR1/CDR2/CDR2.5 contribute unweighted V-loop distance. Profile weights apply
    only to substitutions, in the symmetric longer-junction frame. Native units
    stay common: 6*clipped_cost*position_weight, gap2400*length_difference,
    V-offset200*loop_distance. Exact exclusion punctures identical full encoded
    junctions even for different V calls; duplicate references count separately.
    ``ndn`` uses explicit half-open boundaries: weight10 in V/J and100 in NDN;
    None means uniform100, retaining unavailable boundary cases. Aligned positions
    use the larger query/reference weight. Length-only profiles do not supply boundaries.
    Unknown/ambiguous/mixed-locus inputs fail as in distance_matrix. This supplies
    geometry, not a background population, calibrated P-value or paired count.
    ``v_loop_weight=0`` removes only V-loop cost, retaining identical V eligibility,
    reference/control identity units, CDR3 geometry and exact puncturing. No fitted weight.
    ``same_v_gene=True`` instead conditions on exact resolved V gene identity;
    it preserves the requested V-loop weight and requires thresholds below the native group barrier.

    With ``linear_mass=True``, return ``(counts, masses)`` in one native traversal;
    mass is the integer sum ``max(0, threshold-distance)``. Divide by a positive
    threshold to obtain the sum of triangular kernel weights. Counts include
    boundary ties whose kernel weight is zero. This mode requires the matching
    seqtree native reducer; it never substitutes a dense matrix or another score.
    """
    v_loop_distances = _v_region_context(v_regions, v_loop_distances, species, alphabet, matrix)
    import inspect
    if 'group_distances' not in inspect.signature(gapblock.count_batch).parameters:
        raise RuntimeError('total-distance counts require seqtree gapblock group_distances support')
    if linear_mass and 'linear_mass' not in inspect.signature(gapblock.count_batch).parameters:
        raise RuntimeError('linear mass requires seqtree gapblock linear_mass support')
    if shared_thresholds and ('shared_thresholds' not in inspect.signature(gapblock.count_batch).parameters or
            'shared_thresholds' not in (gapblock._gapblock_count_batch.__doc__ or '')):
        raise RuntimeError('shared cutoffs require Python/native seqtree gapblock shared_thresholds support')
    queries,references,options=_total_distance_options(queries,references,query_v,reference_v,
        species,threads,ctrim,position_weighting,query_boundaries,reference_boundaries,v_loop_weight,same_v_gene,alphabet,stem_depth,stem_weight,matrix,v_loop_distances)
    if same_v_gene and any(t>=2**28 for row in ([thresholds] if shared_thresholds else thresholds) for t in row):
        raise ValueError('conditional V thresholds must be below the native group barrier')
    return gapblock.count_batch(queries,references,thresholds,**options,
        exclude_exact=exclude_exact,**({'linear_mass':True} if linear_mass else {}),
        **({'shared_thresholds':True} if shared_thresholds else {}))


def shortlist_audit_batch(queries,references,query_v,reference_v,*,radius=90,species='human',
                          threads=1,ctrim=2,exclude_exact=False,alphabet='full',stem_depth=0,stem_weight=1.0,matrix='tcrdist',v_loop_distances=None, v_regions=None):
    """Audit the actual CDR3 ball and unchanged V costs in one native reduction.

    Native units, eligibility and full-junction puncturing match total counts.
    No K cap. Older native backends retain a bounded sparse retrieval fallback.
    This audits uniform single-chain retrieval, not annotation confidence.
    """
    v_loop_distances = _v_region_context(v_regions, v_loop_distances, species, alphabet, matrix)
    import time
    import inspect
    streaming=hasattr(gapblock,'group_audit_batch')
    for name in (['group_audit_batch'] if streaming else ['count_batch','ball_batch']):
        required={'position_weights_by_length','group_distances'}
        native=getattr(gapblock,'_gapblock_'+name,None)
        if (not required <= set(inspect.signature(getattr(gapblock,name)).parameters) or
                any(k not in (getattr(native,'__doc__',None) or '') for k in required)):
            raise RuntimeError('shortlist audit requires Python/native positional and grouped audit support')
    if isinstance(radius,bool) or not isinstance(radius,(int,float)) or not np.isfinite(radius) or radius<0:
        raise ValueError('radius must be finite and nonnegative')
    bound=distance_cutoff(radius)
    q,r,options=_total_distance_options(queries,references,query_v,reference_v,
        species,threads,ctrim,'uniform',alphabet=alphabet,stem_depth=stem_depth,stem_weight=stem_weight,matrix=matrix,v_loop_distances=v_loop_distances)
    if streaming:
        start=time.perf_counter()
        before,after,nearest=gapblock.group_audit_batch(q,r,[bound]*len(q),
            exclude_exact=exclude_exact,**options)
        reduced=time.perf_counter()
        return {'cdr3_shortlist_keys':before,'post_v_radius_keys':after,
                'nearest_within_shortlist_distance':[float(v)/200 if v>=0 else None for v in nearest],
                'route':'native_streaming','projected_sparse_payload_bytes':0,
                'native_output_bytes':20*len(q),
                'stage_seconds':{'native_streamed_reduction':reduced-start,
                                 'conversion':time.perf_counter()-reduced}}
    qgroups=np.asarray(options.pop('query_group_ids'),dtype=np.int32)
    rgroups=np.asarray(options.pop('reference_group_ids'),dtype=np.int32)
    offsets=np.asarray(options.pop('group_distances'),dtype=np.int32)
    start=time.perf_counter()
    counts=gapblock.count_batch(q,r,[[bound]]*len(q),exclude_exact=exclude_exact,**options)
    sizes=np.asarray([row[0] for row in counts],dtype=np.int64)
    # ponytail: conservative object estimate; caller guards actual tree RSS; stream natively if too large.
    estimate=int(sizes.sum())*256
    if estimate>24_000_000_000:
        raise MemoryError('CDR3 shortlist projected sparse payload exceeds24GB; no hits truncated')
    counted=time.perf_counter()
    hits=gapblock.ball_batch(q,r,[bound]*len(q),exclude_exact=exclude_exact,**options)
    assert [len(row) for row in hits]==sizes.tolist(),'shortlist count/retrieval differ'
    retrieved=time.perf_counter();n=int(sizes.sum())
    qi=np.repeat(np.arange(len(q)),sizes)
    ri=np.fromiter((h.ref_id for row in hits for h in row),dtype=np.int32,count=n)
    costs=np.fromiter((h.score for row in hits for h in row),dtype=np.int64,count=n)
    del hits
    costs+=offsets[qgroups[qi],rgroups[ri]]
    retained=np.bincount(qi[costs<=bound],minlength=len(q))
    nearest=np.full(len(q),np.iinfo(np.int64).max,dtype=np.int64)
    np.minimum.at(nearest,qi,costs)
    return {'cdr3_shortlist_keys':sizes.tolist(),'post_v_radius_keys':retained.tolist(),
            'nearest_within_shortlist_distance':[float(v)/200 if size else None for v,size in zip(nearest,sizes)],
            'route':'native_ball_then_numpy','projected_sparse_payload_bytes':estimate,'native_output_bytes':None,
            'stage_seconds':{'preflight_count':counted-start,'native_retrieval':retrieved-counted,
                             'vectorized_v':time.perf_counter()-retrieved}}


def _total_distance_options(queries,references,query_v,reference_v,species,threads,ctrim,position_weighting,
                            query_boundaries=None,reference_boundaries=None,v_loop_weight=1,same_v_gene=False,alphabet='full',stem_depth=0,stem_weight=1.0,matrix='tcrdist',v_loop_distances=None):
    """Shared full-coordinate geometry for marginal and linked-pair reductions."""
    if position_weighting not in ('uniform','significance','ndn'):
        raise ValueError('position_weighting must be uniform, significance or ndn')
    if position_weighting!='uniform' and ctrim!=2:
        raise ValueError('full positional geometry requires ctrim default2')
    if type(v_loop_weight) is not int or v_loop_weight not in (0,1):
        raise ValueError("v_loop_weight must be integer0 or1")
    queries, references, model, genes = _distance_inputs(
        queries, references, query_v, reference_v, species, threads, ctrim,v_loop_distances=v_loop_distances)
    costs=_prepared_costs(alphabet,matrix,v_loop_distances)
    loop_dist, indices=_v_loop_groups(genes,model,costs,v_loop_distances)
    lengths={len(s) for s in queries+references}
    if position_weighting=='ndn':
        import inspect
        if ('query_position_weights' not in inspect.signature(gapblock.count_batch).parameters or
            'query_position_weights' not in gapblock._gapblock_count_batch.__doc__):
            raise RuntimeError('NDN geometry requires native per-input positional weights')
        def rows(sequences,boundaries):
            if boundaries is None or len(boundaries)!=len(sequences):
                raise ValueError('each sequence requires an explicit boundary or None fallback')
            out=[]
            for sequence,boundary in zip(sequences,boundaries):
                if boundary is None:
                    out.append([100]*len(sequence))
                    continue
                if len(boundary)!=2 or any(isinstance(x,bool) or not isinstance(x,int) for x in boundary):
                    raise ValueError('boundary endpoints must be integers')
                v,j=boundary
                if not 0<=v<=j<=len(sequence):
                    raise ValueError('boundaries must be half-open junction coordinates')
                out.append([10]*v+[100]*(j-v)+[10]*(len(sequence)-j))
            return out
        sequence_weights=dict(query_position_weights=rows(queries,query_boundaries),
                              reference_position_weights=rows(references,reference_boundaries))
        weights={}
    elif position_weighting=='significance':
        from .regions import significance_weights
        weights={length:[max(1,round(100*w)) for w in significance_weights(length)] for length in lengths}
    else:
        weights={length:[0]*3+[100]*(length-3-ctrim)+[0]*ctrim for length in lengths}
    prepared=v_loop_distances is not None and costs is v_loop_distances['costs']
    native_matrix=(v_loop_distances['count_matrix'] if prepared else
                   SubstitutionMatrix.from_similarity((-3*costs).tolist()))
    scale=(v_loop_distances['distance_units'][0] if prepared else distance_units(costs)[0])
    options=dict(
        matrix=native_matrix,
        gap_open=600*scale,gap_extend=600*scale,gap_prior=_full_gap_prior,
        position_weights_by_length=weights,
        query_group_ids=indices[:len(queries)].tolist(),
        reference_group_ids=indices[len(queries):].tolist(),
        group_distances=(200*v_loop_weight*loop_dist).tolist(),threads=threads)
    if same_v_gene and genes:
        # Preserve allele-loop costs within a gene; reject other genes before alignment.
        alleles=np.unique(genes)
        labels=np.array([g.split('*')[0] for g in alleles])
        options['group_distances']=np.where(labels[:,None]==labels[None,:],
            200*v_loop_weight*loop_dist,2**28).tolist()
    stem_penalty([], [], depth=stem_depth, weight=stem_weight, costs=costs)
    if stem_depth and stem_weight and queries and references:
        import inspect
        required = 'feature_offsets'
        if (not hasattr(gapblock, 'PairedFeatureOffsets') or
                required not in inspect.signature(gapblock.count_batch).parameters or
                required not in (gapblock._gapblock_count_batch.__doc__ or '')):
            raise RuntimeError('stem geometry requires Python/native factored feature offsets')
        if any(len(s) < 2*stem_depth for s in queries+references):
            raise ValueError('stem depth requires junctions of length >= twice the depth')
        index = {a:i for i,a in enumerate(amino_acids())}
        def ends(sequences):
            return [[index[a] for i in range(stem_depth) for a in (s[i],s[-1-i])]
                    for s in sequences]
        options['feature_offsets'] = gapblock.PairedFeatureOffsets(
            ends(queries), ends(references), costs.tolist(), weight=stem_weight, scale=200)
    if position_weighting=='ndn':options.update(sequence_weights)
    return queries,references,options


def _paired_count_geometry(query_alpha,query_beta,reference_alpha,reference_beta,
        query_v_alpha,query_v_beta,reference_v_alpha,reference_v_beta,*,
        species,threads,ctrim,position_weighting,v_loop_weight,same_v_gene,alphabet,
        stem_depth,stem_weight,matrix,v_loop_distances,
        query_boundaries_alpha=None,query_boundaries_beta=None,
        reference_boundaries_alpha=None,reference_boundaries_beta=None):
    """Prepare the identical two-axis geometry for linked native reducers."""
    axes=[];options={};shared=None
    if v_loop_distances is not None and v_loop_distances.get('species') != species:
        raise ValueError('precomputed V loops require matching species')
    model=load_v_loops(species) if v_loop_distances is None else v_loop_distances['model']
    keys=('position_weights_by_length','query_group_ids','reference_group_ids','group_distances')
    if stem_depth and stem_weight:keys+=('feature_offsets',)
    if position_weighting=='ndn':keys+=('query_position_weights','reference_position_weights')
    boundaries={'alpha':(query_boundaries_alpha,reference_boundaries_alpha),
                'beta':(query_boundaries_beta,reference_boundaries_beta)}
    for suffix,q,r,qv,rv in [('alpha',query_alpha,reference_alpha,query_v_alpha,reference_v_alpha),
                            ('beta',query_beta,reference_beta,query_v_beta,reference_v_beta)]:
        locus='TRA' if suffix=='alpha' else 'TRB'
        qv=_resolved_v_calls(qv,model,v_loop_distances,locus=locus)
        rv=_resolved_v_calls(rv,model,v_loop_distances,locus=locus)
        if any(v is None for v in qv+rv):
            raise ValueError('linked pair V calls must resolve to their alpha/beta locus')
        q,r,kw=_total_distance_options(q,r,qv,rv,species,threads,ctrim,position_weighting,*boundaries[suffix],v_loop_weight,same_v_gene,alphabet,stem_depth,stem_weight,matrix,v_loop_distances)
        axes.append((q,r))
        options.update({key+'_'+suffix:kw.pop(key, None) for key in keys})
        shared=kw
    if len(axes[0][0])!=len(axes[1][0]) or len(axes[0][1])!=len(axes[1][1]):
        raise ValueError('each linked pair requires both junctions')
    return axes,shared,options


def paired_total_distance_mass_batch(query_alpha,query_beta,reference_alpha,reference_beta,
        query_v_alpha,query_v_beta,reference_v_alpha,reference_v_beta,thresholds,*,
        species='human',threads=1,ctrim=2,exclude_exact=False,position_weighting='uniform',
        query_boundaries_alpha=None,query_boundaries_beta=None,
        reference_boundaries_alpha=None,reference_boundaries_beta=None,v_loop_weight=1,same_v_gene=False,alphabet='full',stem_depth=0,stem_weight=1.0,matrix='tcrdist',shared_thresholds=False,v_loop_distances=None, v_regions=None):
    """Counts and triangular integer mass for linked alpha+beta total distances.

    Thresholds and masses use native units200 times ordinary TCRdist distance.
    Exact exclusion requires both full junctions equal, irrespective of V calls.
    Uniform full-coordinate trim weights reproduce paired_distance_matrix; no
    dense matrix, independent marginal matching or paired control model is added.
    ``same_v_gene`` requires matching alpha and beta V genes in that same pair,
    preserving the requested V-loop costs; it does not infer joint support from chain marginals.
    """
    v_loop_distances = _v_region_context(v_regions, v_loop_distances, species, alphabet, matrix)
    if not hasattr(gapblock,'paired_sum_count_batch'):
        raise RuntimeError('linked mass requires seqtree paired_sum_count_batch support')
    if shared_thresholds:
        import inspect
        if ('shared_thresholds' not in inspect.signature(gapblock.paired_sum_count_batch).parameters or
                'shared_thresholds' not in (gapblock._gapblock_paired_sum_count_batch.__doc__ or '')):
            raise RuntimeError('shared cutoffs require Python/native seqtree paired_sum_count_batch shared_thresholds support')
    axes,shared,options=_paired_count_geometry(query_alpha,query_beta,reference_alpha,reference_beta,
        query_v_alpha,query_v_beta,reference_v_alpha,reference_v_beta,
        species=species,threads=threads,ctrim=ctrim,position_weighting=position_weighting,
        v_loop_weight=v_loop_weight,same_v_gene=same_v_gene,alphabet=alphabet,stem_depth=stem_depth,
        stem_weight=stem_weight,matrix=matrix,v_loop_distances=v_loop_distances,
        query_boundaries_alpha=query_boundaries_alpha,query_boundaries_beta=query_boundaries_beta,
        reference_boundaries_alpha=reference_boundaries_alpha,reference_boundaries_beta=reference_boundaries_beta)
    if same_v_gene and any(t>=2**28 for row in ([thresholds] if shared_thresholds else thresholds) for t in row):
        raise ValueError('conditional V thresholds must be below the native group barrier')
    return gapblock.paired_sum_count_batch(axes[0][0],axes[1][0],axes[0][1],axes[1][1],thresholds,
        **shared,**options,exclude_exact=exclude_exact,
        **({'shared_thresholds':True} if shared_thresholds else {}))


GEOMETRY_PRESETS=('flank','framev','framemid','framej','freegap','jw','stem4j','jlen','vjlen','norm','jnorm','sw','nw')
J_GEOMETRY_PRESETS=('jw','stem4j','jlen','vjlen','jnorm')


def require_geometry_capabilities(preset, *, paired=False):
    """Check both wrapper and native fields before reading the requested cohort."""
    import inspect
    if preset is None:
        return
    if preset not in GEOMETRY_PRESETS:
        raise ValueError('unknown geometry preset')
    method=getattr(gapblock,'grouped_paired_sum_count_batch' if paired else 'grouped_count_batch',None)
    native=getattr(gapblock,'_GroupedGeometry',None)
    if method is None or native is None:
        raise RuntimeError('geometry extras require updated Python/native grouped support')
    parameters=inspect.signature(method).parameters
    required=[]
    if preset in ('sw','nw'):
        required=[('alignment_mode','alignment_mode')]
    elif preset in ('framev','framemid','framej'):
        required=[('score_quantum_alpha' if paired else 'score_quantum','score_quantum')]
        if paired:required += [('gap_prior_alpha','gap_prior'),('gap_prior_beta','gap_prior')]
    elif preset in J_GEOMETRY_PRESETS:
        required=[('additional_feature_offsets'+suffix,'additional_features')
            for suffix in (('_alpha','_beta') if paired else ('',))]
    if any(parameter not in parameters or not hasattr(native,field) for parameter,field in required):
        raise RuntimeError('geometry extra requires matching updated Python/native '+preset+' support')


def _grouped_preset_options(preset,queries,references,axes,options,context,paired):
    """Fixed recorded diagnostic recipes; alignment and offsets stay factored."""
    costs=context['costs'];scale=int(costs.max())
    if preset in ('sw','nw'):
        options.update(matrix=SubstitutionMatrix.blosum62(),gap_open=max(1,round(scale/2)),
            gap_extend=max(1,round(scale/2)),alignment_mode='local' if preset=='sw' else 'global',
            alignment_scale=600,trim_left=3,trim_right=3)
        options.pop('gap_prior',None)
        for suffix in ('_alpha','_beta') if paired else ('',):
            options.pop('position_weights_by_length'+suffix,None)
    elif preset=='freegap':
        options['gap_prior']=lambda i,d,m:0 if not d or 3<=i<=m-d-3 else gapblock.UNREACHABLE
    for lane,(q,r,qv,rv) in enumerate(axes):
        suffix=('_alpha','_beta')[lane] if paired else ''
        if preset in ('framev','framemid','framej'):
            loci={context['model'][v]['locus'] for v in qv+rv}
            if len(loci)!=1:
                raise ValueError('frame geometry requires one locus')
            column={'TRA':(0,1,3),'TRB':(1,3,6)}[next(iter(loci))][('framev','framemid','framej').index(preset)]
            # Preserve recorded reversed arguments, scale and trimmed coordinates.
            prior=gapblock.frame_prior(column,scale)
            options['gap_prior'+suffix]=lambda i,d,m,prior=prior:(300*prior(i-3,d,m-6)
                if not d or 3<=i<=m-d-3 else gapblock.UNREACHABLE)
            options['score_quantum'+suffix]=600
        if preset=='flank':
            index={a:i for i,a in enumerate(amino_acids())}
            def ends(sequences):
                return [[index[a] for a in s[:3]+s[-3:] for _ in range(2)] for s in sequences]
            options['feature_offsets'+suffix]=gapblock.PairedFeatureOffsets(
                ends(q),ends(r),costs.tolist(),weight=18,scale=200)
        if preset not in J_GEOMETRY_PRESETS:
            continue
        jcol='j'+str(lane)
        if jcol not in queries.columns or jcol not in references.columns:
            raise ValueError('J geometry requires '+jcol+' on both grouped axes')
        js=[('' if value is None else value).split('*')[0]
            for frame in (queries,references) for value in frame[jcol].to_list()]
        terms=[(js,False)]
        if preset in ('jlen','vjlen'):terms.append(([len(s) for s in q+r],True))
        if preset=='vjlen':terms.append(([v.split('*')[0] for v in qv+rv],False))
        descriptors=[]
        for values,length_term in terms:
            keys=list(dict.fromkeys(values));lookup={value:i for i,value in enumerate(keys)}
            ids=[[lookup[value]]*2 for value in values]
            distance=(scale*np.abs(np.asarray(keys)[:,None]-np.asarray(keys)[None,:]) if length_term
                else scale*(np.asarray(keys)[:,None]!=np.asarray(keys)[None,:]))
            descriptors.append(gapblock.PairedFeatureOffsets(ids[:len(q)],ids[len(q):],distance.tolist(),weight=1,scale=200))
        options['additional_feature_offsets'+suffix]=descriptors


def grouped_total_distance_count_batch(queries, references, thresholds, *, paired=False,
        species='human', threads=1, ctrim=2, v_loop_weight=1, alphabet='full',
        stem_depth=0, stem_weight=1.0, matrix='tcrdist', v_loop_distances=None,
        exclude_exact=False, chunk_size=64, v_regions=None, position_weighting='uniform', geometry_preset=None, **memberships):
    """One shared reference traversal with bounded label-count chunks.

    Frames use ``seq0/v0`` (TRA for linked pairs) and ``seq1/v1`` (TRB).
    Label memberships are CSR arrays on the deduplicated reference-key axis.
    Geometry, same-gene auxiliaries and exact puncturing retain the established
    count contract; no kernel arithmetic or query-label fitting occurs here.
    """
    if geometry_preset is not None and geometry_preset not in GEOMETRY_PRESETS:
        raise ValueError('unknown geometry preset')
    if geometry_preset is not None and (ctrim!=3 or alphabet!='full' or position_weighting!='uniform'):
        raise ValueError('geometry presets require full alphabet, trim3 and uniform weights')
    if geometry_preset is not None:
        matrix='tcrdist' if geometry_preset in ('freegap','framev','framemid','framej') else 'blosum62'
        stem_depth=4 if geometry_preset=='stem4j' else 0
        stem_weight=1.0
        if v_loop_distances is not None and v_loop_distances.get('matrix')!=matrix:
            raise ValueError('geometry preset requires matching precomputed matrix')
    require_geometry_capabilities(geometry_preset,paired=paired)
    v_loop_distances = _v_region_context(v_regions, v_loop_distances, species, alphabet, matrix)
    method='grouped_paired_sum_count_batch' if paired else 'grouped_count_batch'
    if not hasattr(gapblock,method):
        raise RuntimeError('grouped kernel counts require updated Python/native seqtree support')
    if position_weighting not in ('uniform','significance-trimmed'):
        raise ValueError('grouped position_weighting must be uniform or significance-trimmed')
    if position_weighting=='significance-trimmed':
        if ctrim!=3:
            raise ValueError('trimmed significance requires trim3 geometry')
        import inspect
        quantum='score_quantum_alpha' if paired else 'score_quantum'
        if (quantum not in inspect.signature(getattr(gapblock,method)).parameters or
                not hasattr(getattr(gapblock,'_GroupedGeometry',None),'score_quantum')):
            raise RuntimeError('trimmed significance requires Python/native alignment quantization support')
    if 'query_withheld_peptide_ids' in memberships:
        import inspect
        if 'query_withheld_peptide_ids' not in inspect.signature(getattr(gapblock,method)).parameters:
            raise RuntimeError('peptide-withheld grouped counts require updated Python/native seqtree support')
        native=getattr(gapblock,'_gapblock_grouped_count_batch',None)
        if not all(name in (getattr(native,'__doc__',None) or '') for name in
                ('query_withheld_peptide_ids','reference_peptide_offsets','reference_peptide_ids')):
            raise RuntimeError('peptide-withheld grouped counts require updated native seqtree support')
    context=(precompute_v_loop_distances(species=species,alphabet=alphabet,matrix=matrix)
             if v_loop_distances is None else v_loop_distances)
    axes=[];gene_axes=[]
    for i in range(2 if paired else 1):
        qv=_resolved_v_calls(queries['v'+str(i)].to_list(),context['model'],context)
        rv=_resolved_v_calls(references['v'+str(i)].to_list(),context['model'],context)
        if any(v is None for v in qv+rv):
            raise ValueError('grouped V calls must resolve before scoring')
        axes.append((queries['seq'+str(i)].to_list(),references['seq'+str(i)].to_list(),qv,rv))
        gene_axes.append([v.split('*')[0] for v in qv+rv])
    import polars as pl
    # Factor the canonical gene tuple once, independent of allele geometry IDs.
    gene_codes=pl.DataFrame({'v'+str(i):values for i,values in enumerate(gene_axes)}).with_columns(
        pl.struct(pl.all()).alias('_gene')).select('_gene').with_row_index('_row')
    codes=gene_codes.select('_gene').unique(maintain_order=True).with_row_index('_gene_id')
    ids=gene_codes.join(codes,on='_gene',how='left',validate='m:1',maintain_order='left')['_gene_id'].to_list()
    if paired:
        axes0,shared,options=_paired_count_geometry(axes[0][0],axes[1][0],axes[0][1],axes[1][1],
            axes[0][2],axes[1][2],axes[0][3],axes[1][3],species=species,threads=threads,
            ctrim=ctrim,position_weighting='uniform',v_loop_weight=v_loop_weight,same_v_gene=False,
            alphabet=alphabet,stem_depth=stem_depth,stem_weight=stem_weight,matrix=matrix,v_loop_distances=context)
        native_axes=(axes0[0][0],axes0[1][0],axes0[0][1],axes0[1][1])
        options={**shared,**options}
    else:
        q,r,options=_total_distance_options(*axes[0],species,threads,ctrim,'uniform',
            v_loop_weight=v_loop_weight,alphabet=alphabet,stem_depth=stem_depth,
            stem_weight=stem_weight,matrix=matrix,v_loop_distances=context)
        native_axes=(q,r)
    if position_weighting=='significance-trimmed':
        from .regions import significance_weights
        for i,(q,r,_,_) in enumerate(axes):
            weights={}
            for length in {len(s) for s in q+r}:
                profile=significance_weights(length)[3:-3]
                mean=sum(profile)/len(profile) if profile else 1.0
                weights[length]=[0]*3+[max(1,round(100*w/mean)) for w in profile]+[0]*3
            suffix=('_alpha','_beta')[i] if paired else ''
            options['position_weights_by_length'+suffix]=weights
            # Native alignment is3 times the prototype raw score: floor to600
            # before adding V/stem offsets equals200*3*(prototype_raw//200).
            options['score_quantum'+suffix]=600
    if geometry_preset is not None:
        _grouped_preset_options(geometry_preset,queries,references,axes,options,context,paired)
    return getattr(gapblock,method)(*native_axes,thresholds,**options,**memberships,
        query_gene_ids=ids[:queries.height],reference_gene_ids=ids[queries.height:],
        exclude_exact=exclude_exact,chunk_size=chunk_size)


def paired_distance_matrix(query_alpha, query_beta, reference_alpha, reference_beta,
                           query_v_alpha, query_v_beta, reference_v_alpha, reference_v_beta,
                           *, species='human', threads=1, ctrim=2, v_loop_weight=1, alphabet='full',
                           stem_depth=0, stem_weight=1.0, matrix='tcrdist',v_loop_distances=None, v_regions=None):
    """Sum alpha/beta distances against the same ordered reference pairs.

    Each position identifies one linked receptor pair; independent marginal nearest
    neighbours must not be substituted for the reference pair distance.
    ``stem_depth`` applies the joint V-J hairpin stem of :func:`stem_penalty` to **both**
    chains. The summed distance carries one stem term per chain. Recorded
    per-chain benchmark views do not establish true linked-pair accuracy. The alpha and beta
    stems are independent hairpins, so summing them is the same composition the rest of
    the paired distance already uses.
    """
    v_loop_distances = _v_region_context(v_regions, v_loop_distances, species, alphabet, matrix)
    parts=list(map(list,(query_alpha,query_beta,reference_alpha,reference_beta,
                        query_v_alpha,query_v_beta,reference_v_alpha,reference_v_beta)))
    qa,qb,ra,rb,qva,qvb,rva,rvb=parts
    if len({len(qa),len(qb),len(qva),len(qvb)})>1 or len({len(ra),len(rb),len(rva),len(rvb)})>1:
        raise ValueError('each linked pair requires alpha and beta junctions and V calls')
    if v_loop_distances is not None and v_loop_distances.get('species') != species:
        raise ValueError('precomputed V loops require matching species')
    model=load_v_loops(species) if v_loop_distances is None else v_loop_distances['model']
    qva,rva=(_resolved_v_calls(v,model,v_loop_distances,locus='TRA') for v in [qva,rva])
    qvb,rvb=(_resolved_v_calls(v,model,v_loop_distances,locus='TRB') for v in [qvb,rvb])
    if any(v is None for v in qva+rva+qvb+rvb):
        raise ValueError('linked pair V calls must resolve uniquely to their declared alpha/beta locus')
    shared=dict(species=species,threads=threads,ctrim=ctrim,v_loop_weight=v_loop_weight,
                alphabet=alphabet,stem_depth=stem_depth,stem_weight=stem_weight,matrix=matrix,v_loop_distances=v_loop_distances)
    out=distance_matrix(qa,ra,qva,rva,**shared)
    out+=distance_matrix(qb,rb,qvb,rvb,**shared)
    return out
