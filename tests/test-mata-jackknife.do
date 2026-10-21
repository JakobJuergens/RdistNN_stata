// test-mata-jackknife.do
// Tests for Mata function rdistnn_jk_cov_from_ranks()
// File under test: src/_rdistnn_mata.ado (rdistnn_jk_cov_from_ranks)
//
// Tags: [smoke] [edge] [bf]

quietly do "src/_rdistnn_mata.ado"
do "tests/test_framework.do"

// ════════════════════════════════════════════════════════════
test_section "rdistnn_jk_cov_from_ranks: [smoke] basic properties"
// ════════════════════════════════════════════════════════════

// [smoke] SE is finite, non-negative
mata {
    rseed(20)
    n = 50; k = 3; s = 8
    _X    = rnormal(n, 2, 0, 1)
    _x0   = (0, 0)
    _rnk  = rdistnn_ranks(_X, _x0, "l2")
    _sig  = rnormal(n, 1, 0, 1)
    _V = rdistnn_jk_cov_from_ranks(_rnk, _sig, k, s, "exact")
    _se = sqrt(max((_V[1,1],0)))
    st_numscalar("_rse",  _se)
    st_numscalar("_rfin", (_se < .) & (_se >= 0))
}
test_assert "[smoke] JK SE is finite and non-negative" scalar(_rfin) == 1

// [smoke] asymp method also returns finite, non-negative SE
mata {
    rseed(21)
    n = 50; k = 3; s = 8
    _X    = rnormal(n, 2, 0, 1)
    _rnk  = rdistnn_ranks(_X, (0, 0), "l2")
    _sig  = rnormal(n, 1, 0, 1)
    _V = rdistnn_jk_cov_from_ranks(_rnk, _sig, k, s, "asymp")
    _se = sqrt(max((_V[1,1],0)))
    st_numscalar("_rfin", (_se < .) & (_se >= 0))
}
test_assert "[smoke] asymp JK SE is finite and non-negative" scalar(_rfin) == 1

// ════════════════════════════════════════════════════════════
test_section "rdistnn_jk_cov_from_ranks: [edge] constant signal"
// ════════════════════════════════════════════════════════════

// [edge] constant signal => all LOO estimates equal theta => SE = 0
mata {
    rseed(22)
    n = 30; k = 2; s = 5
    _X    = rnormal(n, 2, 0, 1)
    _rnk  = rdistnn_ranks(_X, (0, 0), "l2")
    _sig  = J(n, 1, 4.25)    // constant
    _V = rdistnn_jk_cov_from_ranks(_rnk, _sig, k, s, "exact")
    _se = sqrt(max((_V[1,1],0)))
    st_numscalar("_r", _se)
}
test_assert "[edge] constant signal gives SE = 0" scalar(_r) < 1e-10

// ════════════════════════════════════════════════════════════
test_section "rdistnn_jk_cov_from_ranks: [bf] brute-force LOO"
// ════════════════════════════════════════════════════════════

// [bf] Manually compute the delete-1 jackknife and compare.
//   For each i: remove i from the rank vector, shift ranks > r_i down by 1,
//   compute weights on n-1 obs, compute LOO estimate.  Then apply the
//   jackknife variance formula.
mata {
    rseed(23)
    n = 12; k = 2; s = 4
    _X   = rnormal(n, 2, 0, 1)
    _rnk = rdistnn_ranks(_X, (0.1, -0.2), "l2")
    _sig = rnormal(n, 1, 0, 1)

    // fast SE via rdistnn_jk_cov_from_ranks
    _V = rdistnn_jk_cov_from_ranks(_rnk, _sig, k, s, "exact")
    _se_fast = sqrt(max((_V[1,1],0)))

    // brute-force LOO
    _idx = (1::n)
    _loo = J(n, 1, .)
    for (_i = 1; _i <= n; _i++) {
        _kept     = (_idx :!= _i)
        _rnk_m1   = select(_rnk, _kept)
        _rnk_m1   = _rnk_m1 :- (_rnk_m1 :> _rnk[_i])
        _sig_m1   = select(_sig, _kept)
        _w_m1     = rdistnn_weights_from_ranks(_rnk_m1, k, s, "exact")
        _loo[_i]  = sum(_w_m1 :* _sig_m1)
    }
    _se_bf = sqrt((n-1)/n * sum((_loo :- mean(_loo)):^2))

    st_numscalar("_r", abs(_se_fast - _se_bf))
}
test_assert "[bf] fast JK SE matches brute-force LOO" scalar(_r) < 1e-10

// [bf] Same check using asymp weights
mata {
    rseed(24)
    n = 12; k = 2; s = 4
    _X   = rnormal(n, 2, 0, 1)
    _rnk = rdistnn_ranks(_X, (0, 0), "l2")
    _sig = rnormal(n, 1, 0, 1)

    _V = rdistnn_jk_cov_from_ranks(_rnk, _sig, k, s, "asymp")

    _se_fast = sqrt(max((_V[1,1],0)))

    _idx = (1::n)
    _loo = J(n, 1, .)
    for (_i = 1; _i <= n; _i++) {
        _kept     = (_idx :!= _i)
        _rnk_m1   = select(_rnk, _kept)
        _rnk_m1   = _rnk_m1 :- (_rnk_m1 :> _rnk[_i])
        _sig_m1   = select(_sig, _kept)
        _w_m1     = rdistnn_weights_from_ranks(_rnk_m1, k, s, "asymp")
        _loo[_i]  = sum(_w_m1 :* _sig_m1)
    }
    _se_bf = sqrt((n-1)/n * sum((_loo :- mean(_loo)):^2))

    st_numscalar("_r", abs(_se_fast - _se_bf))
}
test_assert "[bf] asymp: fast JK SE matches brute-force LOO" scalar(_r) < 1e-10
