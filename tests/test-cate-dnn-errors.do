// test-cate-dnn-errors.do
// Tests for cate_dnn input validation: all guarded conditions must
// exit with rc=198 (or non-zero for the empty-sample case).
//
// Tags: [errors]

quietly do "src/_rdistnn_mata.ado"
do "tests/test_framework.do"

rdnn_dgp_basic, n(100) seed(5)

// ════════════════════════════════════════════════════════════
test_section "cate_dnn: [errors] required options missing"
// ════════════════════════════════════════════════════════════

// [errors] propensity() is required
cap cate_dnn y, treat(D) x(x1 x2) x0(0 0) k(3) m0(m0) m1(m1)
test_assert "[errors] missing propensity() rejected (rc=198)" _rc == 198

// [errors] at least one of m0() / m1() is required
cap cate_dnn y, treat(D) x(x1 x2) x0(0 0) k(3) propensity(e)
test_assert "[errors] missing both m0() and m1() rejected (rc=198)" _rc == 198

// ════════════════════════════════════════════════════════════
test_section "cate_dnn: [errors] x0 and k/s validation"
// ════════════════════════════════════════════════════════════

// [errors] x0 length mismatch
cap cate_dnn y, treat(D) x(x1 x2) x0(0 0 1) k(3) propensity(e) m0(m0) m1(m1)
test_assert "[errors] x0 length mismatch rejected (rc=198)" _rc == 198

// [errors] k < 1
cap cate_dnn y, treat(D) x(x1 x2) x0(0 0) k(0) propensity(e) m0(m0) m1(m1)
test_assert "[errors] k=0 rejected (rc=198)" _rc == 198

// [errors] s < k
cap cate_dnn y, treat(D) x(x1 x2) x0(0 0) k(6) s(3) ///
    propensity(e) m0(m0) m1(m1)
test_assert "[errors] s < k rejected (rc=198)" _rc == 198

// [errors] s >= n
cap cate_dnn y, treat(D) x(x1 x2) x0(0 0) k(3) s(100) ///
    propensity(e) m0(m0) m1(m1)
test_assert "[errors] s >= n rejected (rc=198)" _rc == 198

// ════════════════════════════════════════════════════════════
test_section "cate_dnn: [errors] distance and method strings"
// ════════════════════════════════════════════════════════════

// [errors] invalid distance string
cap cate_dnn y, treat(D) x(x1 x2) x0(0 0) k(3) ///
    propensity(e) m0(m0) m1(m1) distance(chebyshev)
test_assert "[errors] invalid distance rejected (rc=198)" _rc == 198

// [errors] invalid method string
cap cate_dnn y, treat(D) x(x1 x2) x0(0 0) k(3) ///
    propensity(e) m0(m0) m1(m1) method(bootstrap)
test_assert "[errors] invalid method rejected (rc=198)" _rc == 198

// ════════════════════════════════════════════════════════════
test_section "cate_dnn: [errors] empty sample"
// ════════════════════════════════════════════════════════════

// [errors] no observations after if restriction
cap cate_dnn y if x1 > 1e6, treat(D) x(x1 x2) x0(0 0) k(3) ///
    propensity(e) m0(m0) m1(m1)
test_assert "[errors] empty sample rejected (non-zero rc)" _rc != 0
