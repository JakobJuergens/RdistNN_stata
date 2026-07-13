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

// [errors] both m0() and m1() are required
cap cate_dnn y, treat(D) x(x1 x2) x0(0 0) k(3) propensity(e)
test_assert "[errors] missing both m0() and m1() rejected (rc=198)" _rc == 198

cap cate_dnn y, treat(D) x(x1 x2) x0(0 0) k(3) propensity(e) m0(m0)
test_assert "[errors] missing m1() rejected (rc=198)" _rc == 198

cap cate_dnn y, treat(D) x(x1 x2) x0(0 0) k(3) propensity(e) m1(m1)
test_assert "[errors] missing m0() rejected (rc=198)" _rc == 198

// [errors] outcome and nuisance interface are required without pseudooutcome()
cap cate_dnn, x(x1 x2) x0(0 0)
test_assert "[errors] missing score variable rejected" _rc != 0

cap cate_dnn y, x(x1 x2) x0(0 0) propensity(e) m0(m0) m1(m1)
test_assert "[errors] missing treat() rejected (rc=198)" _rc == 198

cap cate_dnn y, treat(D) x(x1 x2) x0(0 0) propensity(e) mhat(m0) m0(m0) m1(m1)
test_assert "[errors] Robinson and AIPW nuisances rejected together (rc=198)" _rc == 198

cap cate_dnn y, treat(D) x(x1 x2) x0(0 0) propensity(e) m0(m0)
test_assert "[errors] partial AIPW nuisances rejected (rc=198)" _rc == 198

gen double rho_err = y
cap cate_dnn rho_err, x(x1 x2) x0(0 0) pseudooutcome propensity(e)
test_assert "[errors] nuisance rejected with pseudooutcome (rc=198)" _rc == 198

// [errors] exactly one target interface is required
matrix X0_ok = (0,0 \ 0.5,-0.5)
cap cate_dnn y, treat(D) x(x1 x2) propensity(e) m0(m0) m1(m1)
test_assert "[errors] missing x0()/x0matrix() rejected (rc=198)" _rc == 198

cap cate_dnn y, treat(D) x(x1 x2) x0(0 0) x0matrix(X0_ok) ///
    propensity(e) m0(m0) m1(m1)
test_assert "[errors] x0() and x0matrix() rejected together (rc=198)" _rc == 198

cap cate_dnn y, treat(D) x(x1 x2) x0matrix(X0_missing_name) ///
    propensity(e) m0(m0) m1(m1)
test_assert "[errors] nonexistent x0matrix() rejected (rc=198)" _rc == 198

matrix X0_wrong_cols = (0,0,0 \ 1,1,1)
cap cate_dnn y, treat(D) x(x1 x2) x0matrix(X0_wrong_cols) ///
    propensity(e) m0(m0) m1(m1)
test_assert "[errors] x0matrix() column mismatch rejected (rc=198)" _rc == 198

matrix X0_missing = (0,0 \ .,1)
cap cate_dnn y, treat(D) x(x1 x2) x0matrix(X0_missing) ///
    propensity(e) m0(m0) m1(m1)
test_assert "[errors] missing target value rejected (rc=198)" _rc == 198

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

// [errors] even estimate-only fits require s <= n
cap cate_dnn rho_err, x(x1 x2) x0(0 0) k(1) s(101) pseudooutcome nose
test_assert "[errors] s > n rejected with nose (rc=198)" _rc == 198

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
test_section "cate_dnn: [errors] treatment and propensity validation"
// ════════════════════════════════════════════════════════════

// [errors] treatment must be coded 0/1
preserve
replace D = 2 in 1
cap cate_dnn y, treat(D) x(x1 x2) x0(0 0) k(3) propensity(e) m0(m0) m1(m1)
local rc_bad_treat = _rc
restore
test_assert "[errors] non-binary treat() rejected (rc=198)" `rc_bad_treat' == 198

// [errors] propensity must lie in [0,1]
preserve
replace e = -0.1 in 1
cap cate_dnn y, treat(D) x(x1 x2) x0(0 0) k(3) propensity(e) m0(m0) m1(m1)
local rc_low_e = _rc
restore
test_assert "[errors] negative propensity rejected (rc=198)" `rc_low_e' == 198

preserve
replace e = 1.1 in 1
cap cate_dnn y, treat(D) x(x1 x2) x0(0 0) k(3) propensity(e) m0(m0) m1(m1)
local rc_high_e = _rc
restore
test_assert "[errors] propensity above one rejected (rc=198)" `rc_high_e' == 198

preserve
replace e = . in 1
cap cate_dnn y, treat(D) x(x1 x2) x0(0 0) k(3) propensity(e) m0(m0) m1(m1)
local rc_missing_e = _rc
restore
test_assert "[errors] missing propensity rejected (rc=198)" `rc_missing_e' == 198

preserve
replace e = -1e-12 in 1
cap cate_dnn y, treat(D) x(x1 x2) x0(0 0) k(3) propensity(e) m0(m0) m1(m1)
local rc_just_low_e = _rc
restore
test_assert "[errors] propensity just below zero rejected (rc=198)" `rc_just_low_e' == 198

preserve
replace e = 1 + 1e-12 in 1
cap cate_dnn y, treat(D) x(x1 x2) x0(0 0) k(3) propensity(e) m0(m0) m1(m1)
local rc_just_high_e = _rc
restore
test_assert "[errors] propensity just above one rejected (rc=198)" `rc_just_high_e' == 198

foreach bad_clip in 0 .5 -1e-8 . {
    cap cate_dnn y, treat(D) x(x1 x2) x0(0 0) k(3) propensity(e) ///
        m0(m0) m1(m1) clip(`bad_clip')
    test_assert "[errors] invalid clip(`bad_clip') rejected (rc=198)" _rc == 198
}

test_section "cate_dnn: [errors] one-observation Mahalanobis sample"
preserve
keep in 1
cap cate_dnn rho_err, x(x1 x2) x0(0 0) k(1) s(1) ///
    pseudooutcome nose distance(mahalanobis)
local rc_one_mahal_cate = _rc
cap cate_dnn rho_err, x(x1 x2) x0(0 0) k(1) s(1) pseudooutcome nose
local rc_one_l2_cate = _rc
restore
test_assert "[errors] one-observation cate Mahalanobis rejected" `rc_one_mahal_cate' == 2001
test_assert "[endpoint] one-observation cate L2 estimate-only accepted" `rc_one_l2_cate' == 0

// ════════════════════════════════════════════════════════════
test_section "cate_dnn: [errors] empty sample"
// ════════════════════════════════════════════════════════════

// [errors] no observations after if restriction
cap cate_dnn y if x1 > 1e6, treat(D) x(x1 x2) x0(0 0) k(3) ///
    propensity(e) m0(m0) m1(m1)
test_assert "[errors] empty sample rejected (non-zero rc)" _rc != 0
