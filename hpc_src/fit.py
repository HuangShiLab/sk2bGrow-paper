"""PTR estimation from window rates.

Two estimators, for two situations.

**Sorted regression + RANSAC** (:func:`fit_sorted_ransac`) is the iRep/Pilea
skeleton: sort the window log2 rates, regress value on rank, multiply the slope
by the window count. It needs no coordinates at all, which is why it works on
fragmented references — and it is kept here as the parity path for A/B
benchmarking. Its weakness is structural (report defect D3): the fitted range is
set by the largest and smallest window, so the estimate rides on two extreme
order statistics and needs aggressive outlier rejection to stay put.

**V-shape fitting on real coordinates** (:func:`fit_v_shape`) is what TGT
anchors make possible. Anchor positions are known, so log2 coverage can be
regressed directly on circular distance from the origin::

    log2 mu(x) = a - b1 * min(d, k) - b2 * max(0, d - k),   d = dist(x, ori)

With ``b1 == b2`` this is the plain V that CoPTR-Ref shows to be the maximum
likelihood model. The two-slope form exists for multi-fork replication: at
PTR > 2 overlapping rounds put a genuine kink in the profile (report §8.1, R3),
and a single line through a kinked profile is biased. Which of the two is used is
decided by BIC, not asserted.

``ori`` can be supplied (DoriC / Ori-Finder / dnaA) or searched jointly with the
slope. When it is searched, the returned ``ori_confidence`` says how sharply the
data actually pin it down — a slow grower has little gradient and therefore
little information about where the origin is, and the output should say so rather
than reporting a confident coordinate.

Standard errors are scaled by the reduced chi-square of the fit. The window
standard errors from :mod:`sk2bgrow.ztp` describe counting noise only; anchor
efficiency noise, residual GC structure and profile misspecification all add
scatter on top. Taking the residuals at face value keeps the error bars honest,
and those error bars are what :mod:`sk2bgrow.fusion` weights by.

**Poisson GLM on raw window counts** (:func:`fit_glm`, ``method="glm"``) is the
count-native counterpart of the V-shape WLS fit. It models the per-window count
sum directly, with the log2-rate machinery (ZTP truncation, delta-method SEs)
bypassed: ``total_count ~ Poisson(anchors_per_window * mu(x))`` with a log link,
so ``anchors_per_window`` enters as an exposure offset, not a weight. GC
correction enters the same linear predictor as an offset term (natural-log
units), rather than by editing ``log2_rate``. The origin search minimises the
Poisson deviance instead of the WLS SSE — the two criteria coincide only for
pure Gaussian noise. At this stage the GLM is a single-slope V (``b1 == b2``);
the segmented/kink branch of :func:`fit_v_shape` is not ported.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from .gc_bias import tukey_mask

__all__ = [
    "circular_distance",
    "PtrFit",
    "fit_sorted_ransac",
    "fit_v_shape",
    "fit_glm",
    "find_shared_ori",
    "find_shared_ori_glm",
    "fit_windows",
]

_EPS = 1e-12

#: Minimum windows before the two-slope (multi-fork) model may be considered.
MIN_WINDOWS_FOR_SEGMENTED = 30

#: log2(PTR) above which overlapping replication rounds can put a real kink in
#: the profile. Below PTR ~ 2 a kink has no mechanism, so fitting one is
#: over-fitting (design report section 8.1, risk R3).
MULTIFORK_LOG2_PTR = 1.0


def circular_distance(x: np.ndarray, ori: float, genome_len: float) -> np.ndarray:
    """Distance from ``x`` to ``ori`` around a circular chromosome, in bp.

    Ranges over ``[0, genome_len/2]``; the maximum is the terminus.
    """
    d = np.mod(np.asarray(x, dtype=float) - float(ori), float(genome_len))
    return np.minimum(d, float(genome_len) - d)


@dataclass
class PtrFit:
    """One PTR estimate with its uncertainty and provenance."""

    log2_ptr: float
    se: float
    method: str
    n_windows: int
    ori: float | None = None
    ori_confidence: float = np.nan
    r2: float = np.nan
    #: Reduced chi-square. Around 1 means the window standard errors explain the
    #: scatter; >> 1 means extra noise the count model does not capture.
    reduced_chi2: float = np.nan
    #: True when the two-slope (multi-fork) model beat the plain V on BIC.
    segmented: bool = False
    slopes: tuple[float, ...] = field(default_factory=tuple)
    kink: float = np.nan
    ok: bool = True
    note: str = ""

    @property
    def ptr(self) -> float:
        return float(2.0**self.log2_ptr) if np.isfinite(self.log2_ptr) else np.nan


def fit_sorted_ransac(
    log2_rates: np.ndarray,
    n_iter: int = 100,
    trim: float = 0.05,
    inlier_quantile: float = 0.5,
    seed: int = 0,
) -> PtrFit:
    """iRep/Pilea-style sorted regression with RANSAC, for parity benchmarking.

    Windows are sorted by log2 rate and regressed on rank; the slope times the
    window count is log2(PTR). ``n_iter`` RANSAC rounds each fit two random
    points, keep the closest ``inlier_quantile`` of the data, refit, and the
    median slope over rounds is returned — the same construction Pilea uses.

    ``trim`` drops that fraction from each end before fitting, which is the
    minimal defence against D3.
    """
    y = np.asarray(log2_rates, dtype=float)
    y = np.sort(y[np.isfinite(y)])
    n_all = y.size
    if n_all < 5:
        return PtrFit(np.nan, np.nan, "sorted_ransac", n_all, ok=False, note="fewer than 5 usable windows")
    lo = int(np.floor(n_all * trim))
    hi = n_all - lo
    ys = y[lo:hi]
    n = ys.size
    if n < 4:
        return PtrFit(np.nan, np.nan, "sorted_ransac", n_all, ok=False, note="trimming left too few windows")
    x = np.arange(n, dtype=float)

    rng = np.random.default_rng(seed)
    slopes: list[float] = []
    for _ in range(n_iter):
        i, j = rng.choice(n, size=2, replace=False)
        if x[i] == x[j]:
            continue
        m = (ys[j] - ys[i]) / (x[j] - x[i])
        c = ys[i] - m * x[i]
        resid = np.abs(ys - (m * x + c))
        keep = resid <= np.quantile(resid, inlier_quantile)
        if keep.sum() < 3:
            continue
        slopes.append(float(np.polyfit(x[keep], ys[keep], 1)[0]))
    if not slopes:
        return PtrFit(np.nan, np.nan, "sorted_ransac", n_all, ok=False, note="RANSAC found no consensus")

    slope = float(np.median(slopes))
    # The slope is per rank step; the full sorted range spans the whole genome,
    # so scaling by the number of windows recovers the peak-to-trough ratio.
    log2_ptr = slope * n_all
    # Spread across RANSAC rounds is the natural uncertainty for this estimator;
    # 1.4826 * MAD converts it to a standard-deviation scale.
    mad = float(np.median(np.abs(np.array(slopes) - slope)))
    se = 1.4826 * mad * n_all
    fitted = slope * x + float(np.median(ys - slope * x))
    ss_res = float(np.sum((ys - fitted) ** 2))
    ss_tot = float(np.sum((ys - ys.mean()) ** 2))
    return PtrFit(
        log2_ptr=float(log2_ptr),
        se=float(se) if np.isfinite(se) and se > 0 else np.nan,
        method="sorted_ransac",
        n_windows=n_all,
        r2=1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan,
        slopes=(slope,),
        ok=log2_ptr >= 0,
        note="" if log2_ptr >= 0 else "negative slope: no replication gradient",
    )


def _wls(design: np.ndarray, y: np.ndarray, w: np.ndarray) -> tuple[np.ndarray, np.ndarray, float, float]:
    """Weighted least squares. Returns (beta, cov, sse, reduced chi-square)."""
    sw = np.sqrt(w)
    xw = design * sw[:, None]
    yw = y * sw
    beta, *_ = np.linalg.lstsq(xw, yw, rcond=None)
    resid = y - design @ beta
    chi2 = float(np.sum(w * resid**2))
    dof = max(len(y) - design.shape[1], 1)
    red = chi2 / dof
    xtwx = xw.T @ xw
    try:
        inv = np.linalg.inv(xtwx)
    except np.linalg.LinAlgError:
        inv = np.linalg.pinv(xtwx)
    # Scale by the reduced chi-square: the window errors are counting-noise only,
    # so unscaled errors would be optimistic by exactly the factor the residuals
    # reveal.
    return beta, inv * red, float(np.sum(resid**2)), red


def _fit_at(positions, y, w, ori, genome_len, kink):
    """Fit the (possibly segmented) V at a fixed origin and kink."""
    d = circular_distance(positions, ori, genome_len)
    half = genome_len / 2.0
    if kink is None:
        design = np.column_stack([np.ones_like(d), -d])
        contrast = np.array([0.0, half])
    else:
        u1 = np.minimum(d, kink)
        u2 = np.maximum(d - kink, 0.0)
        design = np.column_stack([np.ones_like(d), -u1, -u2])
        contrast = np.array([0.0, kink, half - kink])
    beta, cov, sse, red = _wls(design, y, w)
    log2_ptr = float(contrast @ beta)
    var = float(contrast @ cov @ contrast)
    return beta, cov, sse, red, log2_ptr, max(var, 0.0)


def fit_v_shape(
    positions: np.ndarray,
    log2_rates: np.ndarray,
    genome_len: float,
    se: np.ndarray | None = None,
    ori: float | None = None,
    n_grid: int = 180,
    allow_segmented: bool = True,
    refine: bool = True,
) -> PtrFit:
    """Fit the ori-ter profile on real coordinates.

    ``positions`` are global anchor/window coordinates and ``log2_rates`` the
    matching log2 window rates. ``se`` supplies inverse-variance weights; without
    it every window is weighted equally.

    When ``ori`` is None it is grid-searched then refined. Candidate origins whose
    fit has an *uphill* slope are rejected outright: that solution is the same
    line read from the terminus, and letting it win would report the ter as the
    ori with a negative PTR.
    """
    x = np.asarray(positions, dtype=float)
    y = np.asarray(log2_rates, dtype=float)
    good = np.isfinite(x) & np.isfinite(y)
    x, y = x[good], y[good]
    if se is not None:
        s = np.asarray(se, dtype=float)[good]
        # A window with an unusable error bar keeps a finite but small weight
        # instead of dominating or vanishing.
        w = np.where(np.isfinite(s) & (s > 0), 1.0 / np.maximum(s, _EPS) ** 2, np.nan)
        w = np.where(np.isfinite(w), w, np.nanmedian(w) if np.isfinite(np.nanmedian(w)) else 1.0)
    else:
        w = np.ones_like(y)

    n = y.size
    if n < 5 or genome_len <= 0:
        return PtrFit(np.nan, np.nan, "v_shape", n, ok=False, note="fewer than 5 usable windows")

    half = genome_len / 2.0
    candidates = [float(ori)] if ori is not None else list(np.linspace(0.0, genome_len, n_grid, endpoint=False))

    # The segmented model is gated, not merely BIC-selected. With ~20 windows,
    # two extra parameters buy a BIC improvement by chance often enough that a
    # plain V gets reported as kinked, and the two forms give different
    # log2(PTR) — which then shows up as enzymes disagreeing. The kink is a
    # multi-fork phenomenon, so it is only offered where multi-fork is
    # physically possible (log2 PTR above ~1) and there are enough windows to
    # resolve it.
    kink_grid: list[float | None] = [None]
    if allow_segmented and n >= MIN_WINDOWS_FOR_SEGMENTED:
        probe_ori = float(ori) if ori is not None else float(candidates[0])
        _, _, _, _, plain_ptr, _ = _fit_at(x, y, w, probe_ori, genome_len, None)
        if ori is None:
            # Without a known origin, probe a coarse sweep for the best plain fit
            # before deciding whether a kink is even on the table.
            best_plain = -np.inf
            for o in candidates[:: max(1, len(candidates) // 24)]:
                b, _, _, _, p_, _ = _fit_at(x, y, w, o, genome_len, None)
                if np.all(b[1:] >= 0) and np.isfinite(p_):
                    best_plain = max(best_plain, p_)
            plain_ptr = best_plain
        if np.isfinite(plain_ptr) and plain_ptr >= MULTIFORK_LOG2_PTR:
            kink_grid += [half * f for f in (0.25, 0.4, 0.55, 0.7)]

    best = None
    profile: list[tuple[float, float]] = []
    for o in candidates:
        best_here = np.inf
        for k in kink_grid:
            beta, cov, sse, red, log2_ptr, var = _fit_at(x, y, w, o, genome_len, k)
            slopes = beta[1:]
            if not np.all(np.isfinite(slopes)) or np.any(slopes < 0):
                continue  # uphill: this is the terminus, not the origin
            n_par = design_params(k)
            bic = n * np.log(max(sse / n, _EPS)) + n_par * np.log(n)
            best_here = min(best_here, sse)
            if best is None or bic < best[0]:
                best = (bic, o, k, beta, cov, sse, red, log2_ptr, var)
        profile.append((o, best_here))

    if best is None:
        return PtrFit(np.nan, np.nan, "v_shape", n, ok=False, note="no downhill origin: profile has no replication gradient")

    _, o, k, beta, cov, sse, red, log2_ptr, var = best

    if refine and ori is None:
        # Local refinement around the winning grid point, at 1/20 of the spacing.
        step = genome_len / n_grid
        fine = np.linspace(o - step, o + step, 41) % genome_len
        for o2 in fine:
            b2, c2, s2, r2_, p2, v2 = _fit_at(x, y, w, o2, genome_len, k)
            if np.any(b2[1:] < 0):
                continue
            if s2 < sse:
                o, beta, cov, sse, red, log2_ptr, var = o2, b2, c2, s2, r2_, p2, v2

    ss_tot = float(np.sum((y - y.mean()) ** 2))
    r2 = 1.0 - sse / ss_tot if ss_tot > 0 else np.nan
    conf = _ori_confidence(profile, genome_len, n) if ori is None else 1.0
    return PtrFit(
        log2_ptr=float(log2_ptr),
        se=float(np.sqrt(var)),
        method="v_shape_segmented" if k is not None else "v_shape",
        n_windows=n,
        ori=float(o),
        ori_confidence=float(conf),
        r2=float(r2),
        reduced_chi2=float(red),
        segmented=k is not None,
        slopes=tuple(float(b) for b in beta[1:]),
        kink=float(k) if k is not None else np.nan,
        ok=log2_ptr >= 0,
        note="" if log2_ptr >= 0 else "non-positive PTR",
    )


def design_params(kink: float | None) -> int:
    """Free parameter count: intercept + slope(s), plus the kink when present."""
    return 2 if kink is None else 4


def _ori_confidence(profile: list[tuple[float, float]], genome_len: float, n_obs: int) -> float:
    """Circular mean resultant length of the posterior over origin position.

    The SSE profile is turned into relative likelihood with the Gaussian
    transform ``exp(-n/2 * log(SSE/SSE_min))``, then summarised by the standard
    circular concentration statistic. 1 means the origin is pinned down; near 0
    means the data are consistent with it being almost anywhere — the
    identifiability failure the report flags for slow growers (risk R4).
    """
    vals = np.array([s for _, s in profile], dtype=float)
    oris = np.array([o for o, _ in profile], dtype=float)
    finite = np.isfinite(vals)
    if finite.sum() < 3:
        return 0.0
    vals, oris = vals[finite], oris[finite]
    smin = vals.min()
    if smin <= 0:
        return 0.0
    wts = np.exp(-(n_obs / 2.0) * np.log(vals / smin))
    if not np.isfinite(wts).any() or wts.sum() <= 0:
        return 0.0
    theta = 2.0 * np.pi * oris / genome_len
    cx = float(np.sum(wts * np.cos(theta)) / wts.sum())
    cy = float(np.sum(wts * np.sin(theta)) / wts.sum())
    return float(min(1.0, np.hypot(cx, cy)))


def find_shared_ori(
    windows: pd.DataFrame,
    genome_len: float,
    rate_column: str = "log2_rate",
) -> tuple[float | None, float]:
    """Estimate one origin from all enzymes pooled.

    The origin is a property of the *chromosome*, not of an enzyme. Searching it
    separately per enzyme wastes power — each search sees a fraction of the
    windows — and, worse, injects between-enzyme variance that has nothing to do
    with the biology: two enzymes landing on origins 200 kb apart report
    different slopes, and the cross-enzyme consistency test reads that as the
    enzymes disagreeing when in fact they agree and the *search* diverged.

    Each enzyme's rates are median-centred before pooling, which removes its
    efficiency offset (a per-enzyme intercept) while leaving the shared gradient
    intact.

    Returns ``(ori, confidence)``, or ``(None, nan)`` when there is not enough to
    fit.
    """
    pos, val = [], []
    for _, grp in windows.groupby("enzyme", sort=True):
        y = grp[rate_column].to_numpy(dtype=float)
        finite = np.isfinite(y)
        if finite.sum() < 3:
            continue
        pos.append(grp["global_mid"].to_numpy(dtype=float)[finite])
        val.append(y[finite] - np.median(y[finite]))
    if not pos:
        return None, np.nan
    fit = fit_v_shape(np.concatenate(pos), np.concatenate(val), genome_len, allow_segmented=False)
    if not np.isfinite(fit.log2_ptr) or fit.ori is None:
        return None, np.nan
    return float(fit.ori), float(fit.ori_confidence)


# --------------------------------------------------------------------------
# Poisson GLM estimator (method="glm")
# --------------------------------------------------------------------------

def _poisson_irls(
    design: np.ndarray,
    y: np.ndarray,
    offset: np.ndarray,
    max_iter: int = 60,
    tol: float = 1e-9,
) -> tuple[np.ndarray, np.ndarray, float, float]:
    """IRLS fit of ``log mu = design @ beta + offset`` under a Poisson likelihood.

    Returns ``(beta, cov, deviance, reduced_chi2)``, mirroring the return shape
    of :func:`_wls` so the two estimators stay symmetric. The covariance is
    scaled by the reduced Pearson chi-square, the same residual-scaling policy
    as the WLS path: Poisson is the counting model, and any extra scatter the
    residuals reveal inflates the error bars rather than being ignored.
    """
    x = np.asarray(design, dtype=float)
    y = np.asarray(y, dtype=float)
    off = np.asarray(offset, dtype=float)
    p = x.shape[1]
    beta = np.zeros(p)
    eta = np.clip(off + x @ beta, -30.0, 30.0)  # start mu at the exposure
    for _ in range(max_iter):
        mu = np.exp(eta)
        grad = x.T @ (y - mu)
        h = (x * mu[:, None]).T @ x
        try:
            step = np.linalg.solve(h, grad)
        except np.linalg.LinAlgError:
            break
        if not np.all(np.isfinite(step)):
            break
        eta_new = np.clip(eta + x @ step, -30.0, 30.0)
        if np.max(np.abs(eta_new - eta)) < tol * (1.0 + np.max(np.abs(eta))):
            beta = beta + step
            eta = eta_new
            break
        beta = beta + step
        eta = eta_new
    mu = np.exp(eta)
    y_pos = y > 0
    dev = float(2.0 * np.sum(np.where(y_pos, y * np.log(np.maximum(y, _EPS) / mu), 0.0) - (y - mu)))
    pearson = float(np.sum((y - mu) ** 2 / np.maximum(mu, _EPS)))
    dof = max(len(y) - p, 1)
    red = pearson / dof
    h = (x * mu[:, None]).T @ x
    try:
        inv = np.linalg.inv(h)
    except np.linalg.LinAlgError:
        inv = np.linalg.pinv(h)
    return beta, inv * red, dev, red


def _fit_glm_at(positions, counts, exposure, offset_extra, ori, genome_len):
    """Fit the single-slope V GLM at a fixed origin. Offset/exposure layout and
    return shape mirror :func:`_fit_at`."""
    d = circular_distance(positions, ori, genome_len)
    design = np.column_stack([np.ones_like(d), -d])
    offset = np.log(np.maximum(exposure, _EPS)) + offset_extra
    beta, cov, dev, red = _poisson_irls(design, counts, offset)
    scale = (genome_len / 2.0) / np.log(2.0)
    log2_ptr = float(beta[1] * scale)
    var = float(cov[1, 1] * scale * scale)
    return beta, cov, dev, red, log2_ptr, max(var, 0.0)


def _glm_offset_extra(gc_offset: np.ndarray | None) -> np.ndarray | float:
    """GC term of the GLM linear predictor, in natural-log units.

    :func:`sk2bgrow.gc_bias.apply_to_windows` corrects the V-fit by subtracting
    the window's mean anchor GC offset from ``log2_rate``; the GLM applies the
    identical correction inside its own offset instead (``offset += gc_term``),
    converted to natural log, so the count response is never edited.
    """
    if gc_offset is None:
        return 0.0
    g = np.asarray(gc_offset, dtype=float)
    return -np.log(2.0) * np.where(np.isfinite(g), g, 0.0)


def _deviance_ori_confidence(profile: list[tuple[float, float]], genome_len: float) -> float:
    """Circular mean resultant length of the deviance profile over the origin.

    Same construction as :func:`_ori_confidence`, with the Gaussian SSE
    likelihood ratio ``exp(-n/2 * log(SSE/SSE_min))`` replaced by the Poisson
    likelihood ratio ``exp(-(dev - dev_min) / 2)``.
    """
    vals = np.array([s for _, s in profile], dtype=float)
    oris = np.array([o for o, _ in profile], dtype=float)
    finite = np.isfinite(vals)
    if finite.sum() < 3:
        return 0.0
    vals, oris = vals[finite], oris[finite]
    dmin = vals.min()
    wts = np.exp(-0.5 * (vals - dmin))
    if not np.isfinite(wts).any() or wts.sum() <= 0:
        return 0.0
    theta = 2.0 * np.pi * oris / genome_len
    cx = float(np.sum(wts * np.cos(theta)) / wts.sum())
    cy = float(np.sum(wts * np.sin(theta)) / wts.sum())
    return float(min(1.0, np.hypot(cx, cy)))


def fit_glm(
    positions: np.ndarray,
    total_count: np.ndarray,
    anchors_per_window: np.ndarray,
    genome_len: float,
    gc_offset: np.ndarray | None = None,
    ori: float | None = None,
    n_grid: int = 180,
    refine: bool = True,
) -> PtrFit:
    """Poisson GLM fit of the ori-ter profile on raw per-window counts.

    The response is the window's ``total_count`` (the Poisson sufficient
    statistic that the log2-rate pipeline throws away) and the exposure offset
    is ``log(anchors_per_window)``, so windows enter through the size of the
    count they carry, not through a delta-method error bar. The design reuses
    the V-shape x-axis — piecewise-linear in ``circular_distance`` — but at this
    stage only the single-slope form (``b1 == b2``) is implemented; the
    segmented/kink branch of :func:`fit_v_shape` is not ported.

    ``gc_offset`` is the window's ``mean_gc_offset`` column: it is added to the
    linear predictor as an offset term (natural-log units), the GLM counterpart
    of :func:`sk2bgrow.gc_bias.apply_to_windows`.

    When ``ori`` is None it is grid-searched then refined by minimising the
    Poisson deviance — the GLM's own criterion, not the WLS SSE that
    :func:`fit_v_shape` searches. Uphill slopes are rejected for the same reason
    as in the WLS path: that solution is the terminus read as the origin.
    """
    x = np.asarray(positions, dtype=float)
    y = np.asarray(total_count, dtype=float)
    apw = np.asarray(anchors_per_window, dtype=float)
    good = np.isfinite(x) & np.isfinite(y) & (y >= 0) & np.isfinite(apw) & (apw > 0)
    x, y, apw = x[good], y[good], apw[good]
    if gc_offset is not None:
        gc = _glm_offset_extra(np.asarray(gc_offset, dtype=float)[good])
    else:
        gc = 0.0

    n = y.size
    if n < 5 or genome_len <= 0:
        return PtrFit(np.nan, np.nan, "glm", n, ok=False, note="fewer than 5 usable windows")

    candidates = [float(ori)] if ori is not None else list(np.linspace(0.0, genome_len, n_grid, endpoint=False))

    best = None
    profile: list[tuple[float, float]] = []
    for o in candidates:
        beta, cov, dev, red, log2_ptr, var = _fit_glm_at(x, y, apw, gc, o, genome_len)
        if not np.isfinite(beta[1]) or beta[1] < 0:
            profile.append((o, np.inf))
            continue  # uphill: this is the terminus, not the origin
        profile.append((o, dev))
        if best is None or dev < best[0]:
            best = (dev, o, beta, cov, red, log2_ptr, var)

    if best is None:
        return PtrFit(np.nan, np.nan, "glm", n, ok=False, note="no downhill origin: profile has no replication gradient")

    dev, o, beta, cov, red, log2_ptr, var = best

    if refine and ori is None:
        # Local refinement around the winning grid point, at 1/20 of the spacing.
        step = genome_len / n_grid
        fine = np.linspace(o - step, o + step, 41) % genome_len
        for o2 in fine:
            b2, c2, d2, r2_, p2, v2 = _fit_glm_at(x, y, apw, gc, o2, genome_len)
            if not np.isfinite(b2[1]) or b2[1] < 0:
                continue
            if d2 < dev:
                o, beta, cov, dev, red, log2_ptr, var = o2, b2, c2, d2, r2_, p2, v2

    # Deviance-based R^2 against the intercept-only (constant-rate) model.
    mu0 = max(float(y.sum()) / n, _EPS)
    null_dev = float(2.0 * np.sum(np.where(y > 0, y * np.log(np.maximum(y, _EPS) / mu0), 0.0) - (y - mu0)))
    r2 = 1.0 - dev / null_dev if null_dev > 0 else np.nan
    conf = _deviance_ori_confidence(profile, genome_len) if ori is None else 1.0
    return PtrFit(
        log2_ptr=float(log2_ptr),
        se=float(np.sqrt(var)),
        method="glm",
        n_windows=n,
        ori=float(o),
        ori_confidence=float(conf),
        r2=float(r2),
        reduced_chi2=float(red),
        segmented=False,
        slopes=(float(beta[1]),),
        ok=log2_ptr >= 0,
        note="" if log2_ptr >= 0 else "non-positive PTR",
    )


def find_shared_ori_glm(
    windows: pd.DataFrame,
    genome_len: float,
    n_grid: int = 180,
    refine: bool = True,
) -> tuple[float | None, float]:
    """Estimate one origin from all enzymes pooled, by Poisson deviance.

    The GLM counterpart of :func:`find_shared_ori`: same pooling rationale (the
    origin is a property of the chromosome, and per-enzyme searches inject
    between-enzyme variance that the consistency test would misread as
    disagreement), but each enzyme keeps its own intercept in the linear
    predictor instead of having its rates median-centred. Efficiency offsets are
    then absorbed by the fit rather than pre-subtracted, and the origin is
    chosen by the deviance of the pooled count model.
    """
    parts = []
    for enzyme, grp in windows.groupby("enzyme", sort=True):
        y = grp["total_count"].to_numpy(dtype=float)
        finite = np.isfinite(y) & (y >= 0)
        if finite.sum() < 3:
            continue
        parts.append(
            (
                grp["global_mid"].to_numpy(dtype=float)[finite],
                y[finite],
                grp["anchors_per_window"].to_numpy(dtype=float)[finite],
                _glm_offset_extra(grp["mean_gc_offset"].to_numpy(dtype=float)[finite])
                if "mean_gc_offset" in grp
                else 0.0,
            )
        )
    if not parts:
        return None, np.nan

    n_enz = len(parts)

    def pooled_fit(o: float) -> tuple[float, float] | None:
        """Pooled GLM at origin ``o``: enzyme intercepts + one shared slope."""
        blocks, offs, d_col = [], [], []
        for j, (xj, yj, aj, gj) in enumerate(parts):
            dj = circular_distance(xj, o, genome_len)
            block = np.zeros((len(xj), n_enz))
            block[:, j] = 1.0
            blocks.append(block)
            d_col.append(dj)
            offs.append(np.log(np.maximum(aj, _EPS)) + gj)
        design = np.column_stack([np.concatenate(blocks, axis=0), -np.concatenate(d_col)])
        y_all = np.concatenate([p[1] for p in parts])
        beta, _, dev, _ = _poisson_irls(design, y_all, np.concatenate(offs))
        slope = beta[-1]
        if not np.isfinite(slope) or slope < 0:
            return None
        return dev, slope

    candidates = list(np.linspace(0.0, genome_len, n_grid, endpoint=False))
    best = None
    profile: list[tuple[float, float]] = []
    for o in candidates:
        res = pooled_fit(o)
        if res is None:
            profile.append((o, np.inf))
            continue
        profile.append((o, res[0]))
        if best is None or res[0] < best[0]:
            best = res + (o,)

    if best is None:
        return None, np.nan

    dev, _, o = best
    if refine:
        step = genome_len / n_grid
        for o2 in np.linspace(o - step, o + step, 41) % genome_len:
            res = pooled_fit(o2)
            if res is not None and res[0] < dev:
                dev, o = res[0], o2

    return float(o), _deviance_ori_confidence(profile, genome_len)


def fit_windows(
    windows: pd.DataFrame,
    manifest,
    method: str = "auto",
    min_windows: int = 5,
    tukey_k: float = 1.5,
    rate_column: str = "log2_rate",
    shared_ori: bool = True,
) -> pd.DataFrame:
    """Fit every (sample, genome, enzyme) group and return one row per fit.

    ``method`` is ``"v_shape"``, ``"sorted"``, ``"glm"`` or ``"auto"``.
    ``"auto"`` uses the coordinate fit when the reference is contiguous enough to
    trust a coordinate (see :attr:`sk2bgrow.io.GenomeInfo.is_contiguous`) and
    falls back to sorted regression otherwise — a fragmented MAG has no reliable
    x-axis until ``sk2bgrow scaffold`` has given it one. ``"glm"`` is the Poisson
    GLM on window ``total_count`` (see :func:`fit_glm`).

    ``shared_ori`` estimates one origin per (sample, genome) across all enzymes
    and fits each enzyme's slope at that fixed coordinate — see
    :func:`find_shared_ori` / :func:`find_shared_ori_glm`. Set it to ``False``
    to let every enzyme search independently, which is a useful diagnostic but a
    worse estimator.
    """
    # Resolve one origin per (sample, genome) before fitting any enzyme.
    ori_by_genome: dict[tuple, tuple[float | None, float]] = {}
    if shared_ori:
        for key, grp in windows.groupby(["sample", "genome_id"], sort=True):
            info = manifest.genome(int(key[1]))
            if info.ori is not None:
                ori_by_genome[key] = (float(info.ori), info.ori_confidence)
            elif method != "sorted" and info.is_contiguous:
                if method == "glm":
                    ori_by_genome[key] = find_shared_ori_glm(grp, float(info.genome_len))
                else:
                    ori_by_genome[key] = find_shared_ori(grp, float(info.genome_len), rate_column)

    rows = []
    for (sample, genome_id, genome, enzyme), grp in windows.groupby(
        ["sample", "genome_id", "genome", "enzyme"], sort=True
    ):
        info = manifest.genome(int(genome_id))
        y = grp[rate_column].to_numpy(dtype=float)
        keep = tukey_mask(y, k=tukey_k)
        grp = grp[keep]
        y = y[keep]
        n = int(len(grp))

        ori, ori_conf = ori_by_genome.get(
            (sample, int(genome_id)),
            (float(info.ori) if info.ori is not None else None, info.ori_confidence),
        )

        use_v = method == "v_shape" or (method == "auto" and info.is_contiguous)
        if n < min_windows:
            fit = PtrFit(np.nan, np.nan, "none", n, ok=False, note=f"only {n} windows after outlier trimming")
        elif method == "glm":
            fit = fit_glm(
                grp["global_mid"].to_numpy(dtype=float),
                grp["total_count"].to_numpy(dtype=float),
                grp["anchors_per_window"].to_numpy(dtype=float),
                genome_len=float(info.genome_len),
                gc_offset=grp["mean_gc_offset"].to_numpy(dtype=float) if "mean_gc_offset" in grp else None,
                ori=ori,
            )
            if ori is not None:
                fit.ori = ori
                fit.ori_confidence = ori_conf
        elif use_v:
            fit = fit_v_shape(
                grp["global_mid"].to_numpy(dtype=float),
                y,
                genome_len=float(info.genome_len),
                se=grp["log2_se"].to_numpy(dtype=float) if "log2_se" in grp else None,
                ori=ori,
            )
            if ori is not None:
                fit.ori = ori
                fit.ori_confidence = ori_conf
        else:
            fit = fit_sorted_ransac(y)
            if method == "auto":
                fit.note = (fit.note + "; " if fit.note else "") + f"fragmented reference ({info.n_contigs} contigs)"

        rows.append(
            {
                "sample": sample,
                "genome_id": int(genome_id),
                "genome": genome,
                "enzyme": enzyme,
                "log2_ptr": fit.log2_ptr,
                "ptr": fit.ptr,
                "se": fit.se,
                "method": fit.method,
                "n_windows": fit.n_windows,
                "n_windows_used": n,
                "ori": fit.ori,
                "ori_confidence": fit.ori_confidence,
                "r2": fit.r2,
                "reduced_chi2": fit.reduced_chi2,
                "segmented": fit.segmented,
                "kink": fit.kink,
                "ok": fit.ok,
                "note": fit.note,
                "mean_rate": float(np.nanmean(grp["rate"].to_numpy())) if n else np.nan,
                "mean_detected_fraction": float(np.nanmean(grp["detected_fraction"].to_numpy())) if n else np.nan,
                "mean_dispersion": float(np.nanmean(grp["dispersion"].to_numpy())) if n else np.nan,
                "n_anchors": int(grp["n_anchors"].sum()) if n else 0,
            }
        )
    return pd.DataFrame(rows)
