#!/usr/bin/env python3
"""Estimate, before staging any data, how many callable genotypes a candidate
cohort can be expected to yield at a given set of Y markers.

The estimate is empirical rather than theoretical. It takes a cohort that has
already been mapped and genotyped with this repository's pipeline, measures how
many reads each marker actually attracted per mapped chrY read, and applies that
rate to candidates whose chrY yield is itself predicted from an external
coverage proxy. Nothing here models capture efficiency, probe design or
mappability from first principles: whatever those contribute is already inside
the measured per-marker rate, which is the point.

Two rates are involved and they are calibrated separately:

  1. proxy -> chrY reads.  A candidate's mapped chrY read count is unknown until
     it is mapped, which is the expense the estimate exists to avoid. The
     candidate table therefore supplies a proxy (for the Furtwaengler 2020
     cohort, the `NRY` column of the All Ancient DNA compilation) and this tool
     fits a single through-origin ratio against samples where both the proxy and
     the true mapped count are known.

  2. chrY reads -> reads at marker.  Measured directly from the calibration
     genotype table as (total depth at marker) / (total chrY reads).

Counts at this depth are small — a marker attracting ten reads across fifteen
libraries is a *good* marker here — so the tool reports an exact Poisson
interval on every rate it fits. A power statement built on these numbers that
does not carry the interval is overstating what it knows.

Callability follows the repository's registered call rules, which need >= 2
independent reads before any site is called in either direction. `--min-depth`
defaults to 2 accordingly; it is exposed because the rule is a registered
choice, not a law, and a document that changes it must say so.

The split-power sweep answers a different question from per-sample callability:
given that some number of individuals turn out callable, what is the chance of
observing *both* alleles in the cohort? That is the quantity a prediction of the
form "this cohort should split" actually rests on, and it depends on an unknown
population frequency, so it is swept rather than assumed.
"""

import argparse
import csv
import math
import sys


def poisson_ci(k, lo=0.025, hi=0.975):
    """Exact (Garwood) Poisson confidence interval on a count."""
    import statistics  # noqa: F401  (kept for parity with repo style)

    def chi2_inv(p, df):
        # Wilson-Hilferty inverse, adequate at the precision this is reported to.
        if df <= 0:
            return 0.0
        z = _norm_inv(p)
        return df * (1 - 2.0 / (9 * df) + z * math.sqrt(2.0 / (9 * df))) ** 3

    low = 0.0 if k == 0 else chi2_inv(lo, 2 * k) / 2.0
    high = chi2_inv(hi, 2 * (k + 1)) / 2.0
    return low, high


def _norm_inv(p):
    """Acklam's inverse normal CDF. Sufficient for reporting intervals."""
    a = [-3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02,
         1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00]
    b = [-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02,
         6.680131188771972e+01, -1.328068155288572e+01]
    c = [-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00,
         -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00]
    d = [7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00,
         3.754408661907416e+00]
    pl, ph = 0.02425, 1 - 0.02425
    if p < pl:
        q = math.sqrt(-2 * math.log(p))
        return (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / \
               ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1)
    if p > ph:
        q = math.sqrt(-2 * math.log(1 - p))
        return -(((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / \
               ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1)
    q = p - 0.5
    r = q * q
    return (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q / \
           (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1)


def p_at_least(lam, k):
    """P(Poisson(lam) >= k)."""
    if k <= 0:
        return 1.0
    cum = 0.0
    term = math.exp(-lam)
    for i in range(k):
        cum += term
        term *= lam / (i + 1)
    return 1.0 - cum


def poisson_binomial(ps):
    """Exact distribution of the number of successes, by convolution."""
    dist = [1.0]
    for p in ps:
        nxt = [0.0] * (len(dist) + 1)
        for i, d in enumerate(dist):
            nxt[i] += d * (1 - p)
            nxt[i + 1] += d * p
        dist = nxt
    return dist


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--calib-genotypes", required=True,
                    help="genotype table from y_genotype_batch.sh (sample/marker/depth)")
    ap.add_argument("--calib-coverage", required=True,
                    help="*_coverage.tsv carrying chrY_reads_mq25 per calibration sample")
    ap.add_argument("--candidates", required=True,
                    help="TSV with columns: sample, proxy, and optionally note")
    ap.add_argument("--proxy-calibration", required=True,
                    help="TSV mapping calibration sample -> same proxy, for the ratio fit")
    ap.add_argument("--markers", required=True,
                    help="comma-separated marker names to model")
    ap.add_argument("--min-depth", type=int, default=2,
                    help="depth required for a call under the registered rules (default 2)")
    ap.add_argument("--split-freqs", default="0.1,0.25,0.5,0.75,0.9",
                    help="derived-allele frequencies to sweep for split power")
    ap.add_argument("--rate-scale", type=float, default=1.0,
                    help="multiply every fitted per-marker rate by this factor. The rates are "
                         "fitted on single-digit read counts, so their Poisson interval spans "
                         "a factor of several and dominates every other uncertainty here. "
                         "Re-running at the interval bounds is how a power statement reports "
                         "that honestly instead of quoting the point estimate alone.")
    ap.add_argument("--out-prefix", required=True)
    a = ap.parse_args()

    markers = [m.strip() for m in a.markers.split(",") if m.strip()]

    # --- calibration: chrY reads per sample, and depth per marker -------------
    chry = {}
    for r in read_tsv(a.calib_coverage):
        chry[r["sample"]] = int(r["chrY_reads_mq25"])

    depth = {m: 0 for m in markers}
    seen_samples = set()
    for r in read_tsv(a.calib_genotypes):
        seen_samples.add(r["sample"])
        if r["marker"] in depth:
            depth[r["marker"]] += int(r["depth"])

    # Only calibration samples present in both tables contribute, so the
    # denominator matches the reads the depths were counted from.
    calib_samples = sorted(seen_samples & set(chry))
    total_chry = sum(chry[s] for s in calib_samples)
    if total_chry == 0:
        sys.exit("no chrY reads in calibration set")

    # --- calibration: proxy -> chrY reads ------------------------------------
    proxy_calib = {}
    for r in read_tsv(a.proxy_calibration):
        try:
            v = float(r["proxy"])
        except (ValueError, KeyError):
            continue
        if v > 0:
            proxy_calib[r["sample"]] = v

    pairs = [(proxy_calib[s], chry[s]) for s in calib_samples if s in proxy_calib]
    # Two estimators, both reported. The through-origin fit is the obvious one
    # and is what gets used when the proxy is trustworthy; it is also the one a
    # single corrupt proxy value can move a long way, because a sample with a
    # near-zero proxy and real reads contributes its whole numerator against
    # almost no denominator. The median of per-sample ratios cannot be moved
    # that way, so it is what the prediction uses.
    #
    # This is not a curation step and no sample is dropped: the outlier stays in
    # the table, visible, and simply stops dominating the fit.
    num = sum(y for _, y in pairs)
    den = sum(x for x, _ in pairs)
    ratio_origin = num / den if den else 0.0
    per_sample = sorted(y / x for x, y in pairs)
    if per_sample:
        n = len(per_sample)
        ratio = (per_sample[n // 2] if n % 2
                 else 0.5 * (per_sample[n // 2 - 1] + per_sample[n // 2]))
    else:
        ratio = 0.0

    with open(f"{a.out_prefix}_calibration.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["kind", "key", "value", "ci_low", "ci_high", "note"])
        w.writerow(["proxy_ratio_used", "chrY_per_proxy", f"{ratio:.4f}",
                    f"{per_sample[0]:.4f}" if per_sample else "",
                    f"{per_sample[-1]:.4f}" if per_sample else "",
                    f"median of {len(pairs)} per-sample ratios; ci cols are min/max, not a CI"])
        w.writerow(["proxy_ratio_origin", "chrY_per_proxy", f"{ratio_origin:.4f}", "", "",
                    "through-origin fit, reported for comparison, not used"])
        for s in calib_samples:
            if s in proxy_calib:
                w.writerow(["proxy_sample", s, f"{chry[s] / proxy_calib[s]:.4f}",
                            "", "", f"chrY={chry[s]} proxy={proxy_calib[s]:.0f}"])
        for m in markers:
            k = depth[m]
            lo, hi = poisson_ci(k)
            w.writerow(["marker_rate", m, f"{k / total_chry:.6e}",
                        f"{lo / total_chry:.6e}", f"{hi / total_chry:.6e}",
                        f"{k} reads over {total_chry} chrY reads, {len(calib_samples)} samples"])

    # --- candidates -----------------------------------------------------------
    cands = read_tsv(a.candidates)
    rows = []
    for r in cands:
        try:
            proxy = float(r["proxy"])
        except (ValueError, KeyError):
            proxy = 0.0
        pred_chry = proxy * ratio
        rec = {"sample": r["sample"], "proxy": proxy, "pred_chrY": pred_chry,
               "note": r.get("note", "")}
        for m in markers:
            lam = pred_chry * depth[m] / total_chry * a.rate_scale
            rec[f"lambda_{m}"] = lam
            rec[f"p_call_{m}"] = p_at_least(lam, a.min_depth)
        rows.append(rec)

    with open(f"{a.out_prefix}_candidates.tsv", "w", newline="") as fh:
        cols = (["sample", "proxy", "pred_chrY"] +
                [f"{p}_{m}" for m in markers for p in ("lambda", "p_call")] + ["note"])
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(cols)
        for rec in rows:
            w.writerow([rec["sample"], f"{rec['proxy']:.0f}", f"{rec['pred_chrY']:.0f}"] +
                       [f"{rec[f'{p}_{m}']:.3f}" for m in markers for p in ("lambda", "p_call")] +
                       [rec["note"]])

    # --- cohort-level power ---------------------------------------------------
    freqs = [float(x) for x in a.split_freqs.split(",")]
    with open(f"{a.out_prefix}_power.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["marker", "expected_callable", "p_zero_callable",
                    "p_ge1_callable", "p_ge2_callable", "freq", "p_split"])
        for m in markers:
            ps = [rec[f"p_call_{m}"] for rec in rows]
            dist = poisson_binomial(ps)
            exp = sum(i * d for i, d in enumerate(dist))
            p0 = dist[0]
            p1 = 1 - dist[0]
            p2 = 1 - dist[0] - (dist[1] if len(dist) > 1 else 0.0)
            for f in freqs:
                # Given k callable individuals, a split needs at least one of
                # each allele. Uniformity is not evidence of anything here, so
                # it is deliberately not credited.
                p_split = sum(d * (1 - f ** i - (1 - f) ** i)
                              for i, d in enumerate(dist) if i >= 2)
                w.writerow([m, f"{exp:.2f}", f"{p0:.3f}", f"{p1:.3f}", f"{p2:.3f}",
                            f"{f:.2f}", f"{p_split:.3f}"])

        # Joint callability: the intermediate-branch signature needs a call at
        # two markers in the *same* individual, which is a strictly harder ask
        # than either marker alone and is the quantity that discriminates a real
        # node from a lineage-private one.
        if len(markers) >= 2:
            joint = []
            for rec in rows:
                p = 1.0
                for m in markers:
                    p *= rec[f"p_call_{m}"]
                joint.append(p)
            dist = poisson_binomial(joint)
            exp = sum(i * d for i, d in enumerate(dist))
            w.writerow(["+".join(markers), f"{exp:.2f}", f"{dist[0]:.3f}",
                        f"{1 - dist[0]:.3f}",
                        f"{1 - dist[0] - (dist[1] if len(dist) > 1 else 0):.3f}",
                        "NA", "NA"])

    print(f"proxy ratio chrY/proxy = {ratio:.4f} (median of {len(pairs)}; "
          f"through-origin would give {ratio_origin:.4f})")
    print(f"calibration base = {total_chry} chrY reads over {len(calib_samples)} samples")
    for m in markers:
        lo, hi = poisson_ci(depth[m])
        print(f"  {m}: {depth[m]} reads  rate={depth[m] / total_chry:.3e} "
              f"(95% CI {lo / total_chry:.3e}-{hi / total_chry:.3e})")
    print(f"wrote {a.out_prefix}_{{calibration,candidates,power}}.tsv")


if __name__ == "__main__":
    main()
