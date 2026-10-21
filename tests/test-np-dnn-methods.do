// test-np-dnn-methods.do
// Tests for np_dnn method and distance options.
//
// Tags: [agree] [distance] [opts]

quietly do "src/_rdistnn_mata.ado"
do "tests/test_framework.do"

rdnn_dgp_np, n(300) seed(77)

// ════════════════════════════════════════════════════════════
test_section "np_dnn: [agree] exact vs asymp for large n"
// ════════════════════════════════════════════════════════════

// [agree] With n=300 the exact and asymp point estimates should agree closely.
np_dnn y, xcols(x1 x2) x0(0 0) k(4) s(15) method(exact)
local mu_ex = _b[mu]
local se_ex = _se[mu]

np_dnn y, xcols(x1 x2) x0(0 0) k(4) s(15) method(asymp)
local mu_ap = _b[mu]
local se_ap = _se[mu]

local mu_diff = abs(`mu_ex' - `mu_ap')
local se_diff = abs(`se_ex' - `se_ap')

// L1 norm between estimates should be < 5% of the SE (approximate)
local mu_tol = 0.05 * `se_ex'
test_assert "[agree] exact vs asymp mu within 5% of SE" `mu_diff' < max(`mu_tol', 0.1)
test_assert "[agree] exact vs asymp SE both positive"   `se_ex' > 0 & `se_ap' > 0

// ════════════════════════════════════════════════════════════
test_section "np_dnn: [distance] mahalanobis"
// ════════════════════════════════════════════════════════════

// [distance] mahalanobis runs without error and returns finite results
np_dnn y, xcols(x1 x2) x0(0 0) k(4) s(15) distance(mahalanobis)
local mu_mh = _b[mu]
local se_mh = _se[mu]
local ne_mh = e(neff)

test_assert "[distance] mahal: mu finite"   (`mu_mh' < .)
test_assert "[distance] mahal: se >= 0"      `se_mh' >= 0
test_assert "[distance] mahal: neff > 1"     `ne_mh' > 1

// [distance] mahal e(distance) field is set correctly
local dist_mh "`e(distance)'"
test_assert "[distance] mahal: e(distance) = mahalanobis" "`dist_mh'" == "mahalanobis"

// ════════════════════════════════════════════════════════════
test_section "np_dnn: [opts] method=asymp with mahalanobis distance"
// ════════════════════════════════════════════════════════════

// [opts] combined method=asymp and distance=mahalanobis
np_dnn y, xcols(x1 x2) x0(0 0) k(3) s(10) method(asymp) distance(mahalanobis)
local mu_am = _b[mu]
local se_am = _se[mu]

test_assert "[opts] asymp+mahal: mu finite"   (`mu_am' < .)
test_assert "[opts] asymp+mahal: se >= 0"      `se_am' >= 0

// ════════════════════════════════════════════════════════════
test_section "np_dnn: [agree] reproducibility (two identical calls)"
// ════════════════════════════════════════════════════════════

// [agree] Calling np_dnn twice on the same data with the same arguments
//   must produce identical results (no hidden randomness).
np_dnn y, xcols(x1 x2) x0(0.5 -0.5) k(5) s(20)
local mu1 = _b[mu]
local se1 = _se[mu]

np_dnn y, xcols(x1 x2) x0(0.5 -0.5) k(5) s(20)
local mu2 = _b[mu]
local se2 = _se[mu]

test_assert "[agree] mu is reproducible across identical calls" `mu1' == `mu2'
test_assert "[agree] se is reproducible across identical calls" `se1' == `se2'
