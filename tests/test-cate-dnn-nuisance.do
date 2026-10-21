// test-cate-dnn-nuisance.do
// Tests for cate_dnn nuisance-handling: full AIPW nuisances,
// propensity clipping, oracle recovery, and method/distance options.
//
// Tags: [nuisance] [double-robust] [opts]

quietly do "src/_rdistnn_mata.ado"
do "tests/test_framework.do"

rdnn_dgp_basic, n(200) seed(42)

// ════════════════════════════════════════════════════════════
test_section "cate_dnn: [nuisance] required full AIPW inputs"
// ════════════════════════════════════════════════════════════

// [nuisance] Full AIPW (both m0 and m1 provided)
cate_dnn y, tcol(D) xcols(x1 x2) x0(0 0) k(5) ecol(e) m0col(m0) m1col(m1)
local tau_full = _b[tau]
local se_full  = _se[tau]
test_assert "[nuisance] full AIPW: finite tau"   (`tau_full' < .)
test_assert "[nuisance] full AIPW: se >= 0"       `se_full'  >= 0

// [nuisance] Partial outcome nuisances are not supported.
cap cate_dnn y, tcol(D) xcols(x1 x2) x0(0 0) k(5) ecol(e) m0col(m0)
test_assert "[nuisance] missing m1col() rejected (rc=198)" _rc == 198

cap cate_dnn y, tcol(D) xcols(x1 x2) x0(0 0) k(5) ecol(e) m1col(m1)
test_assert "[nuisance] missing m0col() rejected (rc=198)" _rc == 198

// ════════════════════════════════════════════════════════════
test_section "cate_dnn: [double-robust] oracle nuisances"
// ════════════════════════════════════════════════════════════

// [double-robust] With oracle nuisances and a known simple DGP,
//   the AIPW phi simplifies and the estimate is close to truth.
//   DGP: m0 = 0.5*x2, m1 = sin(x1)+0.5*x2, CATE = sin(x1).
//   At x0=(0,0), true CATE = sin(0) = 0.
rdnn_dgp_basic, n(500) seed(123)
cate_dnn y, tcol(D) xcols(x1 x2) x0(0 0) k(5) ecol(e) m0col(m0) m1col(m1)
local tau_oracle = _b[tau]
local se_oracle  = _se[tau]

// True CATE at (0,0) is sin(0) = 0; estimate should be within 2 SEs
test_assert "[double-robust] oracle: tau within 3 SEs of 0" ///
    abs(`tau_oracle') < 3 * `se_oracle' + 0.3

// ════════════════════════════════════════════════════════════
test_section "cate_dnn: [nuisance] propensity clipping at 1e-8"
// ════════════════════════════════════════════════════════════

// [nuisance] Even with propensity very close to 0 or 1, the command
//   should run without error (residual denominators clip at 1e-8).
rdnn_dgp_basic, n(100) seed(7)
replace e = 1e-9 in 1      // force a near-zero propensity
replace e = 1 - 1e-9 in 2  // force a near-one propensity

cap cate_dnn y, tcol(D) xcols(x1 x2) x0(0 0) k(3) ecol(e) m0col(m0) m1col(m1)
test_assert "[nuisance] propensity clipping: command runs without error" _rc == 0

test_section "cate_dnn: [nuisance] endpoint parity and configurable clipping"

rdnn_dgp_basic, n(100) seed(71)
replace e = 0 in 1
replace e = 1 in 2
replace e = 1e-12 in 3
replace e = 1 - 1e-12 in 4
replace e = .35 in 5
gen double e_used = min(max(e, 1e-6), 1 - 1e-6)
gen double rho_aipw_clip = D*(y-m1)/e_used - (1-D)*(y-m0)/(1-e_used) + m1-m0

quietly cate_dnn y, tcol(D) xcols(x1 x2) x0(0 0) k(3) s(12) ///
    ecol(e) m0col(m0) m1col(m1) clip(1e-6)
local aipw_b = _b[tau]
local aipw_se = _se[tau]
local stored_clip = e(clip)
local propensity_version "`e(propensity_version)'"
quietly cate_dnn rho_aipw_clip, xcols(x1 x2) x0(0 0) k(3) s(12) pseudooutcome
test_assert "[nuisance] AIPW endpoints use configured clipped propensity" ///
    abs(`aipw_b' - _b[tau]) < 1e-12 & abs(`aipw_se' - _se[tau]) < 1e-12
test_assert "[nuisance] clip metadata records estimator value" abs(`stored_clip' - 1e-6) < 1e-15
test_assert "[nuisance] propensity metadata identifies clipped state" ///
    "`propensity_version'" == "clipped"

// ════════════════════════════════════════════════════════════
test_section "cate_dnn: [opts] method and distance options"
// ════════════════════════════════════════════════════════════

rdnn_dgp_basic, n(200) seed(42)

// [opts] method=asymp runs
cate_dnn y, tcol(D) xcols(x1 x2) x0(0 0) k(4) s(15) method(asymp) ///
    ecol(e) m0col(m0) m1col(m1)
local tau_ap = _b[tau]
test_assert "[opts] method=asymp: tau finite" (`tau_ap' < .)

// [opts] distance=mahalanobis runs
cate_dnn y, tcol(D) xcols(x1 x2) x0(0 0) k(4) s(15) distance(mahalanobis) ///
    ecol(e) m0col(m0) m1col(m1)
local tau_mh = _b[tau]
local se_mh  = _se[tau]
test_assert "[opts] distance=mahal: tau finite"  (`tau_mh' < .)
test_assert "[opts] distance=mahal: se >= 0"      `se_mh'  >= 0
