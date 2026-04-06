// test-mata-weights-asymp.do
// Tests for Mata function rdistnn_weights_asymp()
// File under test: src/_rdistnn_mata.ado (rdistnn_weights_asymp)
//
// Tags: [smoke] [edge] [approx] [inv]

quietly do "src/_rdistnn_mata.ado"
do "tests/test_framework.do"

// ════════════════════════════════════════════════════════════
test_section "rdistnn_weights_asymp: [smoke] basic properties"
// ════════════════════════════════════════════════════════════

// [smoke] weights non-negative and sum to 1
mata {
    rseed(10)
    n = 60; k = 4; s = 12
    _rnk = invorder(order(runiform(n, 1), 1))
    _w   = rdistnn_weights_asymp(_rnk, k, s)
    st_numscalar("_rneg", min(_w))
    st_numscalar("_rsum", abs(sum(_w) - 1))
    st_numscalar("_rlen", rows(_w))
}
test_assert "[smoke] asymp weights non-negative"     scalar(_rneg) >= 0
test_assert "[smoke] asymp weights sum to 1"         scalar(_rsum) < 1e-12
test_assert "[smoke] asymp weight vector length = n" scalar(_rlen) == 60

// ════════════════════════════════════════════════════════════
test_section "rdistnn_weights_asymp: [edge] structural cases"
// ════════════════════════════════════════════════════════════

// [edge] s = 1, k = 1: each draw picks one obs uniformly => w_rank = 1/n
//   The Beta(1, 1) = Uniform, so increments are exactly 1/n.
mata {
    n = 31; k = 1; s = 1
    rseed(11)
    _rnk = invorder(order(runiform(n, 1), 1))
    _w   = rdistnn_weights_asymp(_rnk, k, s)
    _wr  = J(n, 1, 0)
    for (i = 1; i <= n; i++) _wr[_rnk[i]] = _w[i]
    st_numscalar("_r", mreldif(_wr, J(n, 1, 1/n)))
}
test_assert "[edge] s=1, k=1: uniform weights 1/n" scalar(_r) < 1e-12

// [edge] k = 1: w_rank should be monotone non-increasing in rank
//   (Beta(1,s) CDF increments on equal bins decrease with rank index)
mata {
    n = 200; k = 1; s = 15
    _rnk = (1::n)   // identity: w[i] = w_rank[i]
    _w   = rdistnn_weights_asymp(_rnk, k, s)
    _diffs = _w[2::n] - _w[1::n-1]   // differences between consecutive w_rank values
    st_numscalar("_r", max(_diffs))   // should be <= 0 (non-increasing)
}
test_assert "[edge] k=1 asymp weights are monotone non-increasing by rank" scalar(_r) <= 1e-12

// ════════════════════════════════════════════════════════════
test_section "rdistnn_weights_asymp: [approx] closeness to exact"
// ════════════════════════════════════════════════════════════

// [approx] For moderate n, asymp and exact L1 distance should be small.
//   Matches the R test: expect l1 < 0.02.
mata {
    n = 120; k = 3; s = 10
    _rnk    = (1::n)
    _w_ex   = rdistnn_weights_exact(_rnk, k, s)
    _w_ap   = rdistnn_weights_asymp(_rnk, k, s)
    _l1     = sum(abs(_w_ex - _w_ap))
    st_numscalar("_r", _l1)
}
test_assert "[approx] asymp vs exact L1 < 0.02 for n=120" scalar(_r) < 0.02

// ════════════════════════════════════════════════════════════
test_section "rdistnn_weights_asymp: [inv] permutation invariance"
// ════════════════════════════════════════════════════════════

// [inv] w_rank depends only on (n,k,s), not on the permutation
mata {
    n = 80; k = 4; s = 12
    rseed(33)

    _r1 = (1::n)
    _r2 = (n::1)
    _r3 = invorder(order(runiform(n, 1), 1))

    _w1 = rdistnn_weights_asymp(_r1, k, s)
    _w2 = rdistnn_weights_asymp(_r2, k, s)
    _w3 = rdistnn_weights_asymp(_r3, k, s)

    _wr1 = J(n, 1, 0); _wr2 = J(n, 1, 0); _wr3 = J(n, 1, 0)
    for (_i = 1; _i <= n; _i++) {
        _wr1[_r1[_i]] = _w1[_i]
        _wr2[_r2[_i]] = _w2[_i]
        _wr3[_r3[_i]] = _w3[_i]
    }

    st_numscalar("_r12", max(abs(_wr2 - _wr1)))
    st_numscalar("_r13", max(abs(_wr3 - _wr1)))
}
test_assert "[inv] asymp w_rank: identity vs reversed"  scalar(_r12) == 0
test_assert "[inv] asymp w_rank: identity vs random"    scalar(_r13) == 0
