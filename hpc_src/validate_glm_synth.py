"""Ad-hoc validation of the Poisson GLM estimator (C8, validation b).

Constructs window count tables with a known origin and known PTR
(single-slope V in natural-log space + Poisson noise, like the prototype),
then checks:

1. fit_glm recovers log2(PTR) within error bars when the origin is supplied.
2. The deviance-based origin search recovers a *shifted* origin (one that has
   nothing to do with the E. coli ori).
3. find_shared_ori_glm recovers the origin when several enzymes with
   different efficiencies are pooled.
4. A window-level GC bias biases the uncorrected fit, and passing
   mean_gc_offset as a GLM offset term removes the bias.
"""

import numpy as np
import pandas as pd

from sk2bgrow.fit import circular_distance, fit_glm, find_shared_ori_glm

GENOME = 4_641_652


def make_counts(true_ori, log2_ptr, n_windows=60, apw=80, intercept=1.2, seed=0):
    rng = np.random.default_rng(seed)
    x = (np.arange(n_windows) + 0.5) / n_windows * GENOME
    d = circular_distance(x, true_ori, GENOME)
    b1 = log2_ptr * np.log(2.0) / (GENOME / 2.0)
    mu_anchor = np.exp(intercept - b1 * d)
    y = rng.poisson(apw * mu_anchor).astype(float)
    return x, y


print("== 1/2. known PTR, shifted origin, ori searched ==")
TRUE_ORI = 1_234_567  # arbitrary offset origin
TRUE_LOG2 = 0.9
x, y = make_counts(TRUE_ORI, TRUE_LOG2, seed=1)
f = fit_glm(x, y, np.full_like(x, 80.0), GENOME)
err_ptr = f.log2_ptr - TRUE_LOG2
err_ori = abs(f.ori - TRUE_ORI)
print(f"true log2_ptr={TRUE_LOG2:.3f}  fitted={f.log2_ptr:.4f} +- {f.se:.4f}  err={err_ptr:+.4f}")
print(f"true ori={TRUE_ORI}  fitted ori={f.ori:.0f}  |err|={err_ori:.0f} bp  conf={f.ori_confidence:.3f}")
print(f"r2={f.r2:.4f} reduced_chi2={f.reduced_chi2:.4f} slopes={f.slopes} method={f.method} ok={f.ok}")
assert f.ok and abs(err_ptr) < max(3 * f.se, 0.02), "PTR recovery failed"
assert err_ori < 60_000, "origin search failed"
assert f.ori_confidence > 0.9, "origin confidence failed"

# Same data, origin supplied: estimate must agree with the searched one.
f2 = fit_glm(x, y, np.full_like(x, 80.0), GENOME, ori=TRUE_ORI)
print(f"ori supplied: log2_ptr={f2.log2_ptr:.4f} +- {f2.se:.4f} (agrees with search: "
      f"{abs(f2.log2_ptr - f.log2_ptr) < 0.02})")
assert abs(f2.log2_ptr - f.log2_ptr) < 0.02

# Consistency: 10x the depth must pin the origin much more tightly.
_, y_deep = make_counts(TRUE_ORI, TRUE_LOG2, apw=800, seed=1)
f_deep = fit_glm(x, y_deep, np.full_like(x, 800.0), GENOME)
err_ori_deep = abs(f_deep.ori - TRUE_ORI)
print(f"10x depth: ori |err|={err_ori_deep:.0f} bp (vs {err_ori:.0f} at 1x), "
      f"log2_ptr err={f_deep.log2_ptr - TRUE_LOG2:+.5f}")
assert err_ori_deep < err_ori / 2, "origin error did not shrink with depth"

print("\n== 3. pooled shared-origin search across 4 enzymes ==")
rows = []
for j, eff in enumerate([1.0, 0.6, 1.7, 0.3]):
    rng = np.random.default_rng(10 + j)
    xj = (np.arange(60) + 0.5) / 60 * GENOME + j * 3000
    dj = circular_distance(xj, TRUE_ORI, GENOME)
    b1 = TRUE_LOG2 * np.log(2.0) / (GENOME / 2.0)
    mu = eff * np.exp(1.2 - b1 * dj)
    yj = rng.poisson(80 * mu).astype(float)
    rows.append(pd.DataFrame({
        "enzyme": f"E{j}", "global_mid": xj, "total_count": yj,
        "anchors_per_window": 80, "mean_gc_offset": np.nan,
    }))
pooled = pd.concat(rows)
o, conf = find_shared_ori_glm(pooled, GENOME)
err_pooled = abs(o - TRUE_ORI)
print(f"pooled ori={o:.0f}  |err|={err_pooled:.0f} bp  conf={conf:.3f}")
# ~15 kb is one sigma here (profile curvature gives sigma ~= 1/sqrt(c) with
# c ~= 8.6e-9/bp^2); allow 3 sigma.
assert err_pooled < 50_000, "pooled origin search failed"

print("\n== 4. GC bias enters as an offset term ==")
rng = np.random.default_rng(7)
x, y = make_counts(TRUE_ORI, TRUE_LOG2, seed=2)
# Window-level GC effect in log2 units (as gc_bias.apply_to_windows would
# subtract it); bake the same factor into the counts.
gc_log2 = 0.6 * np.sin(2 * np.pi * x / GENOME * 3 + 0.5)
y_gc = rng.poisson(np.asarray(y) * 2.0 ** (-gc_log2)).astype(float)
f_raw = fit_glm(x, y_gc, np.full_like(x, 80.0), GENOME, ori=TRUE_ORI)
f_fix = fit_glm(x, y_gc, np.full_like(x, 80.0), GENOME, gc_offset=gc_log2, ori=TRUE_ORI)
print(f"uncorrected: log2_ptr={f_raw.log2_ptr:.4f} (err {f_raw.log2_ptr - TRUE_LOG2:+.4f})")
print(f"offset-corrected: log2_ptr={f_fix.log2_ptr:.4f} +- {f_fix.se:.4f} (err {f_fix.log2_ptr - TRUE_LOG2:+.4f})")
# The offset cancels the baked-in GC factor exactly, so what remains is Poisson
# noise: require the correction to shrink the bias several-fold and the
# residual to sit inside 3 standard errors.
assert abs(f_raw.log2_ptr - TRUE_LOG2) > 3 * abs(f_fix.log2_ptr - TRUE_LOG2), "GC offset had no effect"
assert abs(f_fix.log2_ptr - TRUE_LOG2) < max(3 * f_fix.se, 0.05), "GC offset correction failed"

print("\nALL SYNTHETIC CHECKS PASSED")
