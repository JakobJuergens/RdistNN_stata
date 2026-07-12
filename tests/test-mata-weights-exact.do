// test-mata-weights-exact.do
// Tests for Mata function rdistnn_weights_exact()
// File under test: src/_rdistnn_mata.ado (rdistnn_weights_exact)
//
// Tags: [smoke] [edge] [bf] [inv] [errors]

quietly do "src/_rdistnn_mata.ado"
do "tests/test_framework.do"

// ════════════════════════════════════════════════════════════
test_section "rdistnn_weights_exact: [smoke] basic properties"
// ════════════════════════════════════════════════════════════

// [smoke] weights are non-negative and sum to 1
mata {
    rseed(1)
    n = 12; k = 3; s = 5
    _rnk = invorder(order(runiform(n, 1), 1))   // random permutation of 1..n
    _w   = rdistnn_weights_exact(_rnk, k, s)
    st_numscalar("_rneg",  min(_w))              // should be >= 0
    st_numscalar("_rsum",  abs(sum(_w) - 1))     // should be < eps
    st_numscalar("_rlen",  rows(_w))
}
test_assert "[smoke] weights non-negative"        scalar(_rneg) >= 0
test_assert "[smoke] weights sum to 1"            scalar(_rsum) < 1e-12
test_assert "[smoke] weight vector length = n"    scalar(_rlen) == 12

// ════════════════════════════════════════════════════════════
test_section "rdistnn_weights_exact: [edge] closed-form cases"
// ════════════════════════════════════════════════════════════

// [edge] s = n: subset is always the full sample.
//   i-th NN in any full-sample draw is the i-th global NN.
//   => w_rank[j] = 1/k for j <= k, 0 for j > k.
mata {
    n = 7; k = 3; s = n
    _rnk = (1::n)   // identity permutation
    _w   = rdistnn_weights_exact(_rnk, k, s)
    _exp = (J(k, 1, 1/k) \ J(n-k, 1, 0))
    st_numscalar("_r", max(abs(_w - _exp)))
}
test_assert "[edge] s=n gives uniform 1/k on top-k, 0 elsewhere" scalar(_r) == 0

// [edge] k = 1: closed form is C(n-r, s-1) / C(n, s)
mata {
    n = 10; k = 1; s = 4
    _rnk     = (1::n)
    _w       = rdistnn_weights_exact(_rnk, k, s)
    // reference: lncomb defined in the same Mata session
    _logCns  = lncomb(n, s)
    _w_ref   = J(n, 1, 0)
    for (i = 1; i <= n-s+1; i++) {
        _w_ref[i] = exp(lncomb(n-i, s-1) - _logCns)
    }
    st_numscalar("_r", mreldif(_w, _w_ref))
}
test_assert "[edge] k=1 matches C(n-r,s-1)/C(n,s)" scalar(_r) < 1e-12

// [support] Exact top-k support is ranks 1,...,n-s+k for varied n, s, and k.
mata {
    _support_cases = (10,4,1 \ 12,5,3 \ 20,14,2 \ 9,9,3 \ 9,3,3)
    _support_count_err = 0
    _support_edge_err  = 0
    _support_tail_err  = 0
    for (_case_i = 1; _case_i <= rows(_support_cases); _case_i++) {
        n = _support_cases[_case_i, 1]
        s = _support_cases[_case_i, 2]
        k = _support_cases[_case_i, 3]
        _edge = n - s + k
        _w = rdistnn_weights_exact((1::n), k, s)
        _support_count_err = max((_support_count_err, abs(sum(_w :> 0) - _edge)))
        _support_edge_err  = max((_support_edge_err, _w[_edge] <= 0))
        if (_edge < n) {
            _support_tail_err = max((_support_tail_err, max(abs(_w[(_edge+1)..n]))))
        }
    }
    st_numscalar("_support_count_err", _support_count_err)
    st_numscalar("_support_edge_err",  _support_edge_err)
    st_numscalar("_support_tail_err",  _support_tail_err)
}
test_assert "[support] positive-weight count equals n-s+k" scalar(_support_count_err) == 0
test_assert "[support] rank n-s+k has positive weight" scalar(_support_edge_err) == 0
test_assert "[support] ranks above n-s+k have zero weight" scalar(_support_tail_err) == 0

// [support] Holding n and k fixed, larger s gives tighter support and lower ESS.
mata {
    n = 20; k = 2
    _w_small_s = rdistnn_weights_exact((1::n), k, 5)
    _w_large_s = rdistnn_weights_exact((1::n), k, 12)
    st_numscalar("_support_small_s", sum(_w_small_s :> 0))
    st_numscalar("_support_large_s", sum(_w_large_s :> 0))
    st_numscalar("_ess_small_s", 1 / sum(_w_small_s:^2))
    st_numscalar("_ess_large_s", 1 / sum(_w_large_s:^2))
}
test_assert "[support] larger s shrinks positive support" scalar(_support_large_s) < scalar(_support_small_s)
test_assert "[support] smaller s increases effective local sample size" scalar(_ess_small_s) > scalar(_ess_large_s)

// [edge] s = 1, k = 1: each obs is drawn alone; weight = 1/n uniformly
mata {
    n = 9; k = 1; s = 1
    _rnk = invorder(order(runiform(n, 1), 1))
    _w   = rdistnn_weights_exact(_rnk, k, s)
    // w_rank should all equal 1/n
    _wr  = J(n, 1, 0)
    for (i = 1; i <= n; i++) _wr[_rnk[i]] = _w[i]
    st_numscalar("_r", mreldif(_wr, J(n, 1, 1/n)))
}
test_assert "[edge] s=1, k=1 gives uniform weights 1/n" scalar(_r) < 1e-12

// ════════════════════════════════════════════════════════════
test_section "rdistnn_weights_exact: [bf] brute-force enumeration"
// ════════════════════════════════════════════════════════════

// [bf] Enumerate all C(8,3)=56 subsets; verify exact formula (n=8, s=3, k=2).
//   With identity ranks, the 2 smallest-ranked obs in {i1,i2,i3} (i1<i2<i3)
//   are always i1 and i2.  Each gets weight 1/k per subset.
mata {
    n = 8; s = 3; k = 2
    _rnk  = (1::n)
    _w_ex = rdistnn_weights_exact(_rnk, k, s)

    _acc   = J(n, 1, 0)
    _nsubs = 0
    for (_i1 = 1; _i1 <= n-2; _i1++) {
        for (_i2 = _i1+1; _i2 <= n-1; _i2++) {
            for (_i3 = _i2+1; _i3 <= n; _i3++) {
                _nsubs++
                _acc[_i1] = _acc[_i1] + 1/k
                _acc[_i2] = _acc[_i2] + 1/k
            }
        }
    }
    _w_bf = _acc / _nsubs

    st_numscalar("_rdiff", mreldif(_w_ex, _w_bf))
    st_numscalar("_rsum",  abs(sum(_w_bf) - 1))
}
test_assert "[bf] exact matches brute-force enumeration" scalar(_rdiff) < 1e-10
test_assert "[bf] brute-force weights sum to 1"          scalar(_rsum)  < 1e-12

// ════════════════════════════════════════════════════════════
test_section "rdistnn_weights_exact: [inv] permutation invariance"
// ════════════════════════════════════════════════════════════

// [inv] The weight assigned to rank position r depends only on (n,k,s,r),
//   not on which data observation holds that rank.
//   Equivalently, w_rank (weights reindexed by rank position) is the same
//   regardless of the permutation passed as `ranks'.
mata {
    n = 12; k = 3; s = 5
    rseed(99)

    _r1 = (1::n)                                 // identity
    _r2 = (n::1)                                  // reversed
    _r3 = invorder(order(runiform(n, 1), 1))      // random

    _w1 = rdistnn_weights_exact(_r1, k, s)
    _w2 = rdistnn_weights_exact(_r2, k, s)
    _w3 = rdistnn_weights_exact(_r3, k, s)

    // recover w_rank[r] = w[i] where _rN[i] = r
    _wr1 = J(n, 1, 0); _wr2 = J(n, 1, 0); _wr3 = J(n, 1, 0)
    for (_i = 1; _i <= n; _i++) {
        _wr1[_r1[_i]] = _w1[_i]
        _wr2[_r2[_i]] = _w2[_i]
        _wr3[_r3[_i]] = _w3[_i]
    }

    st_numscalar("_r12", max(abs(_wr2 - _wr1)))
    st_numscalar("_r13", max(abs(_wr3 - _wr1)))
}
test_assert "[inv] w_rank: identity vs reversed permutation"  scalar(_r12) == 0
test_assert "[inv] w_rank: identity vs random permutation"    scalar(_r13) == 0

// Also verify the obs-order mapping: w[i] = w_rank[ranks[i]]
mata {
    n = 10; k = 2; s = 4
    rseed(77)
    _r = invorder(order(runiform(n, 1), 1))
    _w = rdistnn_weights_exact(_r, k, s)

    _wr = J(n, 1, 0)
    for (_i = 1; _i <= n; _i++) _wr[_r[_i]] = _w[_i]

    // w[i] must equal w_rank[ r[i] ]
    _check = J(n, 1, 0)
    for (_i = 1; _i <= n; _i++) _check[_i] = (_w[_i] == _wr[_r[_i]])
    st_numscalar("_r", sum(_check))
}
test_assert "[inv] w[i] = w_rank[ranks[i]] for all i" scalar(_r) == 10

// ════════════════════════════════════════════════════════════
test_section "rdistnn_weights_exact: [errors] out-of-range inputs"
// ════════════════════════════════════════════════════════════

// [errors] k > s should trigger an error
capture mata: _w2 = rdistnn_weights_exact((1::6), 4, 3)
test_assert "[errors] k > s triggers error" _rc != 0

// [errors] s > n should trigger an error
capture mata: _w2 = rdistnn_weights_exact((1::6), 1, 7)
test_assert "[errors] s > n triggers error" _rc != 0
