"""Fixed-K global order-statistic tests, conditional on the declared IID null.

K is declared before search. T is the Kth distinct-key score; control counts
include every score <= T, including ties. These tests do not calibrate labels.
"""

from scipy.stats import beta, binom

from .paired import PAIRED_BOUND_DELTA, _joint_punctured_results


def _upper(count, size, error):
    return 1.0 if count == size else float(beta.ppf(1 - error, count + 1, size - count))


def _validate(k, n, counts):
    if any(isinstance(x, bool) or not isinstance(x, int) for x in (k, n)):
        raise ValueError("K and target exposure must be integers")
    if not 1 <= k <= n:
        raise ValueError("fixed K must lie within the original target exposure")
    for count, size in counts:
        if any(isinstance(x, bool) or not isinstance(x, int) for x in (count, size)):
            raise ValueError("control counts and sizes must be integers")
        if size <= 0 or not 0 <= count <= size:
            raise ValueError("control counts must lie within nonempty control sizes")


def order_statistic(k, n, nc, m):
    """Confidence-bound plus binomial CDF for a predeclared Kth-neighbour score.

    Conditional on the target-selected finite threshold, independent IID controls
    give a pointwise CP upper bound. Adding its failure budget bounds the global
    order-statistic P-value under IID target scores from the same background.
    Deduplicating biological reference/control sequences does not establish IID.
    """
    _validate(k, n, [(nc, m)])
    delta = PAIRED_BOUND_DELTA
    upper = _upper(nc, m, delta)
    return {
        "p_upper": upper,
        "E_raw": n * nc / m,
        "E": n * upper,
        "p_global": min(1.0, delta + float(binom.sf(k - 1, n, upper))),
        "finite_control_delta": delta,
    }


def paired_order_statistic(k, n, na, nb, za, zb, ma, mb, exclude_exact=False):
    """Global fixed-K test for max(chain-A penalty, chain-B penalty).

    Independent background chains define the Cartesian ball. Exact-pair exclusion
    removes its exact/exact corner only; N remains the original pair exposure.
    Sum-score thresholds cannot use this construction.
    """
    _validate(k, n, [(na, ma), (nb, mb), (za, ma), (zb, mb)])
    if za > na or zb > nb:
        raise ValueError("exact control counts must lie inside the selected ball")
    delta = PAIRED_BOUND_DELTA
    if exclude_exact:
        result = _joint_punctured_results([k], [na], [nb], [n], ma, mb, [za], [zb])
        return {
            "p_upper": result["p_upper"][0],
            "E_raw": result["E_raw"][0],
            "E": result["E"][0],
            "p_global": result["p_enrichment"][0],
            "finite_control_delta": delta,
        }
    upper = _upper(na, ma, delta / 2) * _upper(nb, mb, delta / 2)
    return {
        "p_upper": upper,
        "E_raw": n * (na / ma) * (nb / mb),
        "E": n * upper,
        "p_global": min(1.0, delta + float(binom.sf(k - 1, n, upper))),
        "finite_control_delta": delta,
    }
