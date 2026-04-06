// test-cate-dnn-core.do
// Tests for cate_dnn: basic behavior, eclass results, and core properties.
//
// Tags: [smoke] [meta] [edge] [pseudo]

quietly do "src/_rdistnn_mata.ado"
do "tests/test_framework.do"

rdnn_dgp_basic, n(200) seed(42)

// ════════════════════════════════════════════════════════════
test_section "cate_dnn: [smoke] basic run and eclass results"
// ════════════════════════════════════════════════════════════

cate_dnn y, treat(D) x(x1 x2) x0(0 0) k(5) propensity(e) m0(m0) m1(m1)

local tau  = _b[tau]
local se   = _se[tau]
local neff = e(neff)
local N    = e(N)

test_assert "[smoke] _b[tau] is finite"     (`tau' < .)
test_assert "[smoke] _se[tau] is finite"    (`se'  < .)
test_assert "[smoke] _se[tau] >= 0"          `se'  >= 0
test_assert "[smoke] e(neff) > 1"            `neff' > 1
test_assert "[smoke] e(neff) <= e(N)"        `neff' <= `N'
test_assert "[smoke] e(N) = 200"             `N'   == 200

// ════════════════════════════════════════════════════════════
test_section "cate_dnn: [meta] eclass fields"
// ════════════════════════════════════════════════════════════

cate_dnn y, treat(D) x(x1 x2) x0(0 0) k(3) s(15) ///
    propensity(e) m0(m0) m1(m1) distance(l2) method(exact)

local cmd  "`e(cmd)'"
local dist "`e(distance)'"
local meth "`e(method)'"
local k_r  = e(k)
local s_r  = e(s)
local kx_r = e(kx)

test_assert "[meta] e(cmd) = cate_dnn"       "``cmd''"  == "cate_dnn"
test_assert "[meta] e(distance) = l2"        "``dist''" == "l2"
test_assert "[meta] e(method) = exact"       "``meth''" == "exact"
test_assert "[meta] e(k) = 3"               `k_r'  == 3
test_assert "[meta] e(s) = 15"              `s_r'  == 15
test_assert "[meta] e(kx) = 2"              `kx_r' == 2

// ════════════════════════════════════════════════════════════
test_section "cate_dnn: [pseudo] tau = weighted mean of AIPW phi"
// ════════════════════════════════════════════════════════════

// [pseudo] Manually compute full AIPW phi and DNN-weighted mean;
//   result must equal _b[tau] exactly (up to floating-point).
rdnn_dgp_basic, n(200) seed(42)

// Compute AIPW pseudo-outcome manually
gen double den  = max(e * (1 - e), 1e-8)
gen double phi  = D * (y - m1) / den  ///
                - (1-D) * (y - m0) / den  ///
                + (m1 - m0)

// Get DNN weights
weights_dnn x1 x2, x0(0 0) k(5) generate(w_cate)

// Manual tau
mata {
    _w   = st_data(., "w_cate")
    _phi = st_data(., "phi")
    st_numscalar("_tau_manual", sum(_w :* _phi))
}
local tau_manual = scalar(_tau_manual)

// Command tau
cate_dnn y, treat(D) x(x1 x2) x0(0 0) k(5) ///
    propensity(e) m0(m0) m1(m1)
local tau_cmd = _b[tau]

test_assert "[pseudo] tau = weighted mean of AIPW phi" abs(`tau_cmd' - `tau_manual') < 1e-8

// ════════════════════════════════════════════════════════════
test_section "cate_dnn: [edge] constant pseudo-outcome => SE = 0"
// ════════════════════════════════════════════════════════════

// [edge] If all pseudo-outcomes are constant, the LOO JK SE is 0.
rdnn_dgp_basic, n(100) seed(11)
replace m0 = y          // set m0 = y so phi = constant (m1-m0) = constant
replace m1 = y + 2.0    // tau_true = 2.0

cate_dnn y, treat(D) x(x1 x2) x0(0 0) k(3) ///
    propensity(e) m0(m0) m1(m1)
local se_const = _se[tau]
test_assert "[edge] constant phi => se = 0" `se_const' < 1e-8
