"""Exact induced-subgraph acquisition under uniform sampling without replacement."""
import numpy as np
import polars as pl
from scipy.special import gammaln


def degree_diversity(degrees):
    """Degree-based Hill summaries and an explicitly labelled Chao count proxy.

    For uniform vertex mass p=1/N and binary self-inclusive Hamming similarity
    Z, ordinariness (Zp)_i=(d_i+1)/N yields similarity-aware Hill orders 0,1,2.
    The separate Chao adaptation treats d+1 as a synthetic count vector; it is
    not an unseen-family estimator under the original record-sampling model.
    """
    from vdjtools.stats import asymptotic_diversity, sample_coverage

    d = np.asarray(degrees)
    if d.ndim != 1 or d.dtype.kind not in 'iu' or (d < 0).any() or (d >= len(d)).any() or int(d.sum()) % 2:
        raise ValueError('nonnegative integer graph degrees below N required')
    if not len(d):
        return dict(degree_diversity_status='empty_graph', degree_chao_status='empty_graph')
    counts = d+1
    z = counts/len(d)
    entropy = float(-np.log(z).mean())
    hills = [float(np.mean(1/z)), float(np.exp(entropy)), float(1/z.mean())]
    return dict(degree_diversity_status='uniform_vertex_binary_similarity',
        graph_similarity_entropy_nats=entropy,
        **{f'graph_hill{q}': value for q, value in enumerate(hills)},
        **{f'graph_hill{q}_enrichment': len(d)/value for q, value in enumerate(hills)},
        degree_mass=int(counts.sum()), degree_singletons=int((counts == 1).sum()),
        degree_doubletons=int((counts == 2).sum()),
        degree_chao_coverage=float(sample_coverage(counts)),
        degree_chao1=float(asymptotic_diversity(counts, q=0, se=False)['estimator'][0]),
        degree_chao_status='degree_mass_proxy')


def graph_acquisition(degrees, sizes):
    """Expected edges and connected sampled vertices of the observed fixed graph.

    ``degrees`` excludes self and represents distinct junction vertices. For
    sample size m from N vertices, a sampled degree-d vertex is orphaned with
    probability choose(N-1-d,m-1)/choose(N-1,m-1). This is interpolation of the
    observed graph, not unseen-family coverage or Chao's taxon-abundance model.
    """
    d, sizes = np.asarray(degrees), np.asarray(sizes)
    if (d.ndim != 1 or sizes.ndim != 1 or not np.issubdtype(d.dtype, np.integer)
            or not np.issubdtype(sizes.dtype, np.integer) or (d < 0).any()
            or (d >= len(d)).any() or int(d.sum()) % 2
            or (sizes < 1).any() or (sizes > len(d)).any()):
        raise ValueError('integer graph degrees and sample sizes within 1..N required')
    N = len(d)
    values, counts = np.unique(d, return_counts=True)
    edges = int(d.sum())//2
    rows = []
    for m in np.unique(sizes):
        # Binomial ratios cancel factorial(m-1); only degree histogram bins
        # enter this reduction, never an N by N distance or adjacency matrix.
        k, a, b = int(m)-1, N-1-values, N-1
        orphan = np.zeros(len(values))
        valid = a >= k
        orphan[valid] = np.exp(np.minimum(0., gammaln(a[valid]+1)-gammaln(a[valid]-k+1)
            -gammaln(b+1)+gammaln(b-k+1)))
        fraction = float(np.dot(orphan, counts)/N)
        rows.append(dict(graph_junctions=N, sample_junctions=int(m),
            expected_graph_edges=edges*int(m)*(int(m)-1)/(N*(N-1)) if N > 1 else 0.,
            expected_orphan_junctions=int(m)*fraction,
            expected_connected_junctions=int(m)*(1-fraction),
            expected_connected_fraction=1-fraction))
    return pl.DataFrame(rows, schema={'graph_junctions': pl.Int64, 'sample_junctions': pl.Int64,
        'expected_graph_edges': pl.Float64, 'expected_orphan_junctions': pl.Float64,
        'expected_connected_junctions': pl.Float64, 'expected_connected_fraction': pl.Float64})
