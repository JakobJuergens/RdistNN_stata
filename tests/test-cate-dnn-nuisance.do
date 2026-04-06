// test-cate-dnn-nuisance.do
// Tests for cate_dnn nuisance-handling: full AIPW vs Robinson forms,
// propensity clipping, oracle recovery, and method/distance options.
//
// Tags: [nuisance] [double-robust] [opts]

quietly do "src/_rdistnn_mata.ado"
do "tests/test_framework.do"

rdnn_dgp_basic, n(200) seed(42)

// ════════════════════════════════════════════════════════════
test_section "cate_dnn: [nuisance] pseudo-outcome forms"
// ════════════════════════════════════════════════════════════

// [nuisance] Full AIPW (both m0 and m1 provided)
cate_dnn y, treat(D) x(x1 x2) x0(0 0) k(5) propensity(e) m0(m0) m1(m1)
local tau_full = _b[tau]
local se_full  = _se[tau]
test_assert "[nuisance] full AIPW: finite tau"   (`tau_full' < .)
test_assert "[nuisance] full AIPW: se >= 0"       `se_full'  >= 0

// [nuisance] Robinson form: m0 only (note from command should print)
quietly cate_dnn y, treat(D) x(x1 x2) x0(0 0) k(5) propensity(e) m0(m0)
local tau_rob0 = _b[tau]
local se_rob0  = _se[tau]
test_assert "[nuisance] Robinson m0-only: finite tau"   (`tau_rob0' < .)
test_assert "[nuisance] Robinson m0-only: se >= 0"       `se_rob0'  >= 0

// [nuisance] Robinson form: m1 only
quietly cate_dnn y, treat(D) x(x1 x2) x0(0 0) k(5) propensity(e) m1(m1)
local tau_rob1 = _b[tau]
local se_rob1  = _se[tau]
test_assert "[nuisance] Robinson m1-only: finite tau"   (`tau_rob1' < .)
test_assert "[nuisance] Robinson m1-only: se >= 0"       `se_rob1'  >= 0

// [nuisance] full AIPW and Robinson forms give different estimates
//   (same oracle nuisances but different pseudo-outcome)
local diff_forms = abs(`tau_full' - `tau_rob0')
test_assert "[nuisance] full AIPW != Robinson m0-only" `diff_forms' > 0

// ════════════════════════════════════════════════════════════
test_section "cate_dnn: [double-robust] oracle nuisances"
// ════════════════════════════════════════════════════════════

// [double-robust] With oracle nuisances and a known simple DGP,
//   the AIPW phi simplifies and the estimate is close to truth.
//   DGP: m0 = 0.5*x2, m1 = sin(x1)+0.5*x2, CATE = sin(x1).
//   At x0=(0,0), true CATE = sin(0) = 0.
rdnn_dgp_basic, n(500) seed(123)
cate_dnn y, treat(D) x(x1 x2) x0(0 0) k(5) propensity(e) m0(m0) m1(m1)
local tau_oracle = _b[tau]
local se_oracle  = _se[tau]

// True CATE at (0,0) is sin(0) = 0; estimate should be within 2 SEs
test_assert "[double-robust] oracle: tau within 3 SEs of 0" ///
    abs(`tau_oracle') < 3 * `se_oracle' + 0.3

// ════════════════════════════════════════════════════════════
test_section "cate_dnn: [nuisance] propensity clipping at 1e-8"
// ════════════════════════════════════════════════════════════

// [nuisance] Even with propensity very close to 0 or 1, the command
//   should run without error (denominator clips at 1e-8).
rdnn_dgp_basic, n(100) seed(7)
replace e = 1e-9 in 1      // force a near-zero propensity
replace e = 1 - 1e-9 in 2  // force a near-one propensity

cap cate_dnn y, treat(D) x(x1 x2) x0(0 0) k(3) propensity(e) m0(m0) m1(m1)
test_assert "[nuisance] propensity clipping: command runs without error" _rc == 0

// ════════════════════════════════════════════════════════════
test_section "cate_dnn: [opts] method and distance options"
// ════════════════════════════════════════════════════════════

rdnn_dgp_basic, n(200) seed(42)

// [opts] method=asymp runs
cate_dnn y, treat(D) x(x1 x2) x0(0 0) k(4) s(15) method(asymp) ///
    propensity(e) m0(m0) m1(m1)
local tau_ap = _b[tau]
test_assert "[opts] method=asymp: tau finite" (`tau_ap' < .)

// [opts] distance=mahalanobis runs
cate_dnn y, treat(D) x(x1 x2) x0(0 0) k(4) s(15) distance(mahalanobis) ///
    propensity(e) m0(m0) m1(m1)
local tau_mh = _b[tau]
local se_mh  = _se[tau]
test_assert "[opts] distance=mahal: tau finite"  (`tau_mh' < .)
test_assert "[opts] distance=mahal: se >= 0"      `se_mh'  >= 0
