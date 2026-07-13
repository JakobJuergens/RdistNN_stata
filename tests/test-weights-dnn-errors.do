// test-weights-dnn-errors.do
// Tests for weights_dnn input validation: every guarded condition should
// exit with rc=198.
//
// Tags: [errors]

quietly do "src/_rdistnn_mata.ado"
do "tests/test_framework.do"

rdnn_dgp_basic, n(100) seed(1)

// ════════════════════════════════════════════════════════════
test_section "weights_dnn: [errors] input validation"
// ════════════════════════════════════════════════════════════

// [errors] k < 1 must be rejected
cap weights_dnn x1 x2, x0(0 0) k(0)
test_assert "[errors] k=0 rejected (rc=198)" _rc == 198

// [errors] x0 length mismatch: 2 covariates but 3 x0 values
cap weights_dnn x1 x2, x0(0 0 1) k(3)
test_assert "[errors] x0 length mismatch rejected (rc=198)" _rc == 198

// [errors] x0 length mismatch: 2 covariates but 1 x0 value
cap weights_dnn x1 x2, x0(0) k(3)
test_assert "[errors] x0 too short rejected (rc=198)" _rc == 198

// [errors] invalid distance string
cap weights_dnn x1 x2, x0(0 0) k(3) distance(euclidean)
test_assert "[errors] invalid distance string rejected (rc=198)" _rc == 198

// [errors] invalid method string
cap weights_dnn x1 x2, x0(0 0) k(3) method(bruteforce)
test_assert "[errors] invalid method string rejected (rc=198)" _rc == 198

// [errors] s < k: explicit s smaller than k
cap weights_dnn x1 x2, x0(0 0) k(5) s(3)
test_assert "[errors] s < k rejected (rc=198)" _rc == 198

// [endpoint] s = n is valid for standalone weight construction
cap weights_dnn x1 x2, x0(0 0) k(3) s(100)
test_assert "[endpoint] s = n accepted" _rc == 0

// [errors] s > n (explicit)
cap weights_dnn x1 x2, x0(0 0) k(3) s(200)
test_assert "[errors] s > n rejected (rc=198)" _rc == 198

// [errors] mahalanobis with single covariate (need >= 2 for cov matrix)
cap weights_dnn x1, x0(0) k(3) distance(mahalanobis)
test_assert "[errors] mahal with kx=1 rejected (rc=198)" _rc == 198

preserve
keep in 1
cap weights_dnn x1 x2, x0(0 0) k(1) s(1) distance(mahalanobis)
local rc_one_mahal_weights = _rc
cap weights_dnn x1 x2, x0(0 0) k(1) s(1) generate(w_one_l2)
local rc_one_l2_weights = _rc
local w_one = w_one_l2[1]
restore
test_assert "[errors] one-observation weights Mahalanobis rejected" `rc_one_mahal_weights' == 2001
test_assert "[endpoint] one-observation exact L2 weights accepted" ///
    `rc_one_l2_weights' == 0 & abs(`w_one' - 1) < 1e-12
