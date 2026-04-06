// test-np-dnn-smoke.do
// Tests for np_dnn: basic behavior, eclass results, metadata, and edge cases.
//
// Tags: [smoke] [meta] [edge]

quietly do "src/_rdistnn_mata.ado"
do "tests/test_framework.do"

rdnn_dgp_np, n(200) seed(7)

// ════════════════════════════════════════════════════════════
test_section "np_dnn: [smoke] basic run and eclass results"
// ════════════════════════════════════════════════════════════

np_dnn y, x(x1 x2) x0(0 0) k(5)

local mu   = _b[mu]
local se   = _se[mu]
local neff = e(neff)
local N    = e(N)

test_assert "[smoke] _b[mu] is finite"          (`mu' < .)
test_assert "[smoke] _se[mu] is finite"         (`se' < .)
test_assert "[smoke] _se[mu] >= 0"              `se' >= 0
test_assert "[smoke] e(neff) > 1"               `neff' > 1
test_assert "[smoke] e(neff) <= e(N)"           `neff' <= `N'
test_assert "[smoke] e(N) = 200"                `N'    == 200

// [smoke] mu is in plausible range: DGP truth at (0,0) is 0 - 0 = 0
//   With n=200, the estimate should be within 1 of the truth.
test_assert "[smoke] mu at (0,0) is plausible" abs(`mu') < 1

// ════════════════════════════════════════════════════════════
test_section "np_dnn: [meta] eclass scalar and string fields"
// ════════════════════════════════════════════════════════════

np_dnn y, x(x1 x2) x0(0 0) k(3) s(15) distance(l2) method(exact)
local cmd  "`e(cmd)'"
local dist "`e(distance)'"
local meth "`e(method)'"
local k_r  = e(k)
local s_r  = e(s)
local kx_r = e(kx)

test_assert "[meta] e(cmd) = np_dnn"         "``cmd''" == "np_dnn"
test_assert "[meta] e(distance) = l2"        "``dist''" == "l2"
test_assert "[meta] e(method) = exact"       "``meth''" == "exact"
test_assert "[meta] e(k) = 3"               `k_r'  == 3
test_assert "[meta] e(s) = 15"              `s_r'  == 15
test_assert "[meta] e(kx) = 2"              `kx_r' == 2

// ════════════════════════════════════════════════════════════
test_section "np_dnn: [edge] constant outcome"
// ════════════════════════════════════════════════════════════

// [edge] When y is constant, _b[mu] = constant and _se[mu] = 0
rdnn_dgp_const, n(100) seed(1) yval(4.25)

np_dnn y, x(x1 x2) x0(0 0) k(3)
local mu_const = _b[mu]
local se_const = _se[mu]

test_assert "[edge] constant y: mu = yval"        abs(`mu_const' - 4.25) < 1e-8
test_assert "[edge] constant y: se = 0"            `se_const'            < 1e-8

// ════════════════════════════════════════════════════════════
test_section "np_dnn: [smoke] estimate = weighted mean"
// ════════════════════════════════════════════════════════════

// [smoke] _b[mu] equals sum(w_i * y_i) where w comes from weights_dnn
rdnn_dgp_np, n(150) seed(99)

weights_dnn x1 x2, x0(0 0) k(4) s(10) generate(w_np)
quietly summarize w_np    // won't interfere with following ereturn
mata {
    _w = st_data(., "w_np")
    _y = st_data(., "y")
    st_numscalar("_mu_manual", sum(_w :* _y))
}
local mu_manual = scalar(_mu_manual)

np_dnn y, x(x1 x2) x0(0 0) k(4) s(10)
local mu_cmd = _b[mu]

test_assert "[smoke] mu = sum(w*y)" abs(`mu_cmd' - `mu_manual') < 1e-8
