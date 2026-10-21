// test-np-dnn-errors.do
// Tests for np_dnn input validation: all guarded conditions must exit rc=198.
//
// Tags: [errors]

quietly do "src/_rdistnn_mata.ado"
do "tests/test_framework.do"

rdnn_dgp_np, n(100) seed(3)

// ════════════════════════════════════════════════════════════
test_section "np_dnn: [errors] input validation"
// ════════════════════════════════════════════════════════════

// [errors] x0 length mismatch: 2 covariates, 3 x0 values
cap np_dnn y, xcols(x1 x2) x0(0 0 1) k(3)
test_assert "[errors] x0 length mismatch rejected (rc=198)" _rc == 198

// [errors] x0 too short: 2 covariates, 1 x0 value
cap np_dnn y, xcols(x1 x2) x0(0) k(3)
test_assert "[errors] x0 too short rejected (rc=198)" _rc == 198

// [errors] k < 1
cap np_dnn y, xcols(x1 x2) x0(0 0) k(0)
test_assert "[errors] k=0 rejected (rc=198)" _rc == 198

// [errors] s < k (explicit s)
cap np_dnn y, xcols(x1 x2) x0(0 0) k(6) s(4)
test_assert "[errors] s < k rejected (rc=198)" _rc == 198

// [errors] s >= n
cap np_dnn y, xcols(x1 x2) x0(0 0) k(3) s(100)
test_assert "[errors] s >= n rejected (rc=198)" _rc == 198

// [errors] invalid distance string
cap np_dnn y, xcols(x1 x2) x0(0 0) k(3) distance(minkowski)
test_assert "[errors] invalid distance rejected (rc=198)" _rc == 198

// [errors] invalid method string
cap np_dnn y, xcols(x1 x2) x0(0 0) k(3) method(mc)
test_assert "[errors] invalid method rejected (rc=198)" _rc == 198

// [errors] no observations after if restriction
cap np_dnn y if x1 > 1e6, xcols(x1 x2) x0(0 0) k(3)
test_assert "[errors] empty sample rejected (non-zero rc)" _rc != 0

preserve
keep in 1
cap np_dnn y, xcols(x1 x2) x0(0 0) k(1) s(1) nose distance(mahalanobis)
local rc_one_mahal_np = _rc
cap np_dnn y, xcols(x1 x2) x0(0 0) k(1) s(1) nose
local rc_one_l2_np = _rc
restore
test_assert "[errors] one-observation NP Mahalanobis rejected" `rc_one_mahal_np' == 2001
test_assert "[endpoint] one-observation NP L2 estimate-only accepted" `rc_one_l2_np' == 0
