// test-plqr-dnn-errors.do
// Tests for plqr_dnn input validation and missing-data handling.
//
// Tags: [plqr] [errors] [missing]

quietly do "src/_rdistnn_mata.ado"
do "tests/test_framework.do"

rdnn_dgp_plqr, n(100) seed(601)

// ════════════════════════════════════════════════════════════
test_section "plqr_dnn: [errors] required options"
// ════════════════════════════════════════════════════════════

cap plqr_dnn y, x(x1 x2) x0(0 0) q(q) md(md1 md2) beta(1 -0.5)
test_assert "[errors] missing d() rejected (rc=198)" _rc == 198

cap plqr_dnn y, x(x1 x2) x0(0 0) d(d1 d2) md(md1 md2) beta(1 -0.5)
test_assert "[errors] missing q() rejected (rc=198)" _rc == 198

cap plqr_dnn y, x(x1 x2) x0(0 0) d(d1 d2) q(q) beta(1 -0.5)
test_assert "[errors] missing md() rejected (rc=198)" _rc == 198

cap plqr_dnn y, x(x1 x2) x0(0 0) d(d1 d2) q(q) md(md1 md2)
test_assert "[errors] missing beta() rejected (rc=198)" _rc == 198

// ════════════════════════════════════════════════════════════
test_section "plqr_dnn: [errors] shape and tuning validation"
// ════════════════════════════════════════════════════════════

cap plqr_dnn y, x(x1 x2) x0(0) d(d1 d2) q(q) md(md1 md2) beta(1 -0.5)
test_assert "[errors] x0 length mismatch rejected (rc=198)" _rc == 198

cap plqr_dnn y, x(x1 x2) x0(0 0) d(d1 d2) q(q) md(md1 md2) beta(1)
test_assert "[errors] beta length mismatch rejected (rc=198)" _rc == 198

cap plqr_dnn y, x(x1 x2) x0(0 0) d(d1 d2) q(q) md(md1) beta(1 -0.5)
test_assert "[errors] d()/md() length mismatch rejected (rc=198)" _rc == 198

cap plqr_dnn y, x(x1 x2) x0(0 0) d(d1 d2) q(q) md(md1 md2) beta(1 -0.5) ///
    tau(0)
test_assert "[errors] tau=0 rejected (rc=198)" _rc == 198

cap plqr_dnn y, x(x1 x2) x0(0 0) d(d1 d2) q(q) md(md1 md2) beta(1 -0.5) ///
    tau(1)
test_assert "[errors] tau=1 rejected (rc=198)" _rc == 198

cap plqr_dnn y, x(x1 x2) x0(0 0) d(d1 d2) q(q) md(md1 md2) beta(1 -0.5) ///
    k(0)
test_assert "[errors] k=0 rejected (rc=198)" _rc == 198

cap plqr_dnn y, x(x1 x2) x0(0 0) d(d1 d2) q(q) md(md1 md2) beta(1 -0.5) ///
    k(5) s(3)
test_assert "[errors] s < k rejected (rc=198)" _rc == 198

cap plqr_dnn y, x(x1 x2) x0(0 0) d(d1 d2) q(q) md(md1 md2) beta(1 -0.5) ///
    k(3) s(100)
test_assert "[errors] s >= n rejected (rc=198)" _rc == 198

// ════════════════════════════════════════════════════════════
test_section "plqr_dnn: [errors] method and distance validation"
// ════════════════════════════════════════════════════════════

cap plqr_dnn y, x(x1 x2) x0(0 0) d(d1 d2) q(q) md(md1 md2) beta(1 -0.5) ///
    distance(chebyshev)
test_assert "[errors] invalid distance rejected (rc=198)" _rc == 198

cap plqr_dnn y, x(x1 x2) x0(0 0) d(d1 d2) q(q) md(md1 md2) beta(1 -0.5) ///
    method(bootstrap)
test_assert "[errors] invalid method rejected (rc=198)" _rc == 198

// ════════════════════════════════════════════════════════════
test_section "plqr_dnn: [missing] markout excludes missing inputs"
// ════════════════════════════════════════════════════════════

rdnn_dgp_plqr, n(100) seed(602)
replace y = . in 1
replace d1 = . in 2
replace md2 = . in 3
replace q = . in 4
plqr_dnn y, x(x1 x2) x0(0 0) d(d1 d2) q(q) md(md1 md2) beta(1 -0.5) ///
    tau(0.5) k(3) s(10)
local N_markout = e(N)
test_assert "[missing] missing y/d/md/q rows are excluded" `N_markout' == 96

replace y = .
cap plqr_dnn y, x(x1 x2) x0(0 0) d(d1 d2) q(q) md(md1 md2) beta(1 -0.5) ///
    tau(0.5) k(3) s(10)
test_assert "[missing] all-missing estimation sample rejected" _rc == 2000

// ════════════════════════════════════════════════════════════
test_section "plqr_dnn: [errors] empty if sample"
// ════════════════════════════════════════════════════════════

rdnn_dgp_plqr, n(100) seed(603)
cap plqr_dnn y if x1 > 1e6, x(x1 x2) x0(0 0) d(d1 d2) q(q) md(md1 md2) ///
    beta(1 -0.5) tau(0.5) k(3) s(10)
test_assert "[errors] empty if sample rejected" _rc == 2000
