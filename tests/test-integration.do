// test-integration.do
// End-to-end pipeline tests: verify that the three public commands are
// mutually consistent and that the DNN estimator recovers known values
// in well-controlled DGPs.
//
// Tags: [e2e] [agree] [consistency]

quietly do "src/_rdistnn_mata.ado"
do "tests/test_framework.do"

// ════════════════════════════════════════════════════════════
test_section "integration: [e2e] weights -> np pipeline consistency"
// ════════════════════════════════════════════════════════════

// [e2e] np_dnn estimate must equal sum(w_i * y_i) where w comes from
//   weights_dnn with the same k, s, x0, method.
rdnn_dgp_np, n(200) seed(10)

weights_dnn x1 x2, x0(0.2 -0.3) k(4) s(12) generate(w_e2e)

mata {
    _w = st_data(., "w_e2e")
    _y = st_data(., "y")
    st_numscalar("_mu_w", sum(_w :* _y))
}
local mu_w = scalar(_mu_w)

np_dnn y, xcols(x1 x2) x0(0.2 -0.3) k(4) s(12)
local mu_cmd = _b[mu]

test_assert "[e2e] np_dnn mu = sum(w*y) from weights_dnn" abs(`mu_cmd' - `mu_w') < 1e-8

// ════════════════════════════════════════════════════════════
test_section "integration: [e2e] cate_dnn with zero treatment effect"
// ════════════════════════════════════════════════════════════

// [e2e] If treatment has no effect (m0 = m1), AIPW phi reduces to
//   D*(y-m0)/e - (1-D)*(y-m0)/(1-e), which has conditional mean 0.
//   With large n and oracle nuisances, |tau_hat| should be small.
clear
set obs 500
set seed 99
gen double x1 = rnormal()
gen double x2 = rnormal()
gen double e  = invlogit(0.5 * x1 - 0.3 * x2)
gen byte   D  = (runiform() < e)
gen double m0 = x1 - 0.5 * x2    // same for D=0 and D=1 (no effect)
gen double m1 = m0                // identical outcome models
gen double y  = m0 + rnormal() * 0.5

cate_dnn y, tcol(D) xcols(x1 x2) x0(0 0) k(5) ecol(e) m0col(m0) m1col(m1)
local tau_null = _b[tau]
local se_null  = _se[tau]

// tau should be close to 0; allow up to 3 SEs
test_assert "[e2e] zero-effect DGP: |tau| < 3*se + 0.2" ///
    abs(`tau_null') < 3 * `se_null' + 0.2

// ════════════════════════════════════════════════════════════
test_section "integration: [agree] all three commands on same data"
// ════════════════════════════════════════════════════════════

// [agree] weights_dnn, np_dnn, and cate_dnn all use the same Mata core.
//   Verify that neff returned by each command agrees.
rdnn_dgp_basic, n(200) seed(42)

weights_dnn x1 x2, x0(0 0) k(5) s(20)
local neff_w = r(neff)

np_dnn y, xcols(x1 x2) x0(0 0) k(5) s(20)
local neff_np = e(neff)

cate_dnn y, tcol(D) xcols(x1 x2) x0(0 0) k(5) s(20) ecol(e) m0col(m0) m1col(m1)
local neff_ct = e(neff)

// All three should return identical neff (same weights, same evaluation point)
test_assert "[agree] weights_dnn neff = np_dnn neff"    abs(`neff_w' - `neff_np') < 1e-8
test_assert "[agree] weights_dnn neff = cate_dnn neff"  abs(`neff_w' - `neff_ct') < 1e-8

// ════════════════════════════════════════════════════════════
test_section "integration: [consistency] ereturn survives replay"
// ════════════════════════════════════════════════════════════

// [consistency] Calling np_dnn with no arguments replays the last results.
rdnn_dgp_np, n(150) seed(20)
np_dnn y, xcols(x1 x2) x0(1 0) k(4)
local mu1 = _b[mu]
local se1 = _se[mu]

np_dnn    // replay
local mu2 = _b[mu]
local se2 = _se[mu]

test_assert "[consistency] replay: mu unchanged" `mu1' == `mu2'
test_assert "[consistency] replay: se unchanged" `se1' == `se2'
