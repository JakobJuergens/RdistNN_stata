// test-mata-jackknife-vector.do
// Tests for Mata function rdistnn_jk_cov_from_ranks()
// File under test: src/_rdistnn_mata.ado
//
// Tags: [jackknife] [matrix] [bf]

quietly do "src/_rdistnn_mata.ado"
do "tests/test_framework.do"

// ════════════════════════════════════════════════════════════
test_section "rdistnn_jk_cov_from_ranks: [compat] scalar signal"
// ════════════════════════════════════════════════════════════

mata {
    rseed(301)
    n = 30; k = 2; s = 5
    _X   = rnormal(n, 2, 0, 1)
    _rnk = rdistnn_ranks(_X, (0.1, -0.2), "l2")
    _sig = rnormal(n, 1, 0, 1)

    _V_exact  = rdistnn_jk_cov_from_ranks(_rnk, _sig, k, s, "exact")
    _se_exact = rdistnn_jk_se_from_ranks(_rnk, _sig, k, s, "exact")
    _V_asymp  = rdistnn_jk_cov_from_ranks(_rnk, _sig, k, s, "asymp")
    _se_asymp = rdistnn_jk_se_from_ranks(_rnk, _sig, k, s, "asymp")

    st_numscalar("_r_exact", abs(_V_exact[1,1] - _se_exact^2))
    st_numscalar("_r_asymp", abs(_V_asymp[1,1] - _se_asymp^2))
}
test_assert "[compat] exact covariance equals scalar SE squared" scalar(_r_exact) < 1e-10
test_assert "[compat] asymp covariance equals scalar SE squared" scalar(_r_asymp) < 1e-10

// ════════════════════════════════════════════════════════════
test_section "rdistnn_jk_cov_from_ranks: [edge] constant matrix signal"
// ════════════════════════════════════════════════════════════

mata {
    rseed(302)
    n = 28; k = 2; s = 6
    _X   = rnormal(n, 2, 0, 1)
    _rnk = rdistnn_ranks(_X, (0, 0), "l2")
    _sig = J(n, 1, 1) * (2.5, -1.25)
    _V   = rdistnn_jk_cov_from_ranks(_rnk, _sig, k, s, "exact")
    st_numscalar("_r_const", max(abs(_V)))
}
test_assert "[edge] constant matrix signal gives zero covariance" scalar(_r_const) < 1e-10

// ════════════════════════════════════════════════════════════
test_section "rdistnn_jk_cov_from_ranks: [bf] brute-force covariance"
// ════════════════════════════════════════════════════════════

mata {
    rseed(303)
    n = 14; q = 3; k = 2; s = 5
    _X   = rnormal(n, 2, 0, 1)
    _rnk = rdistnn_ranks(_X, (0.2, -0.4), "l2")
    _sig = rnormal(n, q, 0, 1)

    _V_fast = rdistnn_jk_cov_from_ranks(_rnk, _sig, k, s, "exact")

    _idx = (1::n)
    _loo = J(n, q, .)
    for (_i = 1; _i <= n; _i++) {
        _kept      = (_idx :!= _i)
        _rnk_m1    = select(_rnk, _kept)
        _rnk_m1    = _rnk_m1 :- (_rnk_m1 :> _rnk[_i])
        _sig_m1    = select(_sig, _kept)
        _w_m1      = rdistnn_weights_from_ranks(_rnk_m1, k, s, "exact")
        _loo[_i,.] = cross(_w_m1, _sig_m1)
    }
    _centered = _loo :- mean(_loo)
    _V_bf = ((n - 1) / n) * cross(_centered, _centered)

    st_numscalar("_r_exact_bf", max(abs(_V_fast - _V_bf)))
}
test_assert "[bf] exact vector covariance matches brute-force LOO" scalar(_r_exact_bf) < 1e-10

mata {
    rseed(304)
    n = 14; q = 3; k = 2; s = 5
    _X   = rnormal(n, 2, 0, 1)
    _rnk = rdistnn_ranks(_X, (0, 0), "l2")
    _sig = rnormal(n, q, 0, 1)

    _V_fast = rdistnn_jk_cov_from_ranks(_rnk, _sig, k, s, "asymp")

    _idx = (1::n)
    _loo = J(n, q, .)
    for (_i = 1; _i <= n; _i++) {
        _kept      = (_idx :!= _i)
        _rnk_m1    = select(_rnk, _kept)
        _rnk_m1    = _rnk_m1 :- (_rnk_m1 :> _rnk[_i])
        _sig_m1    = select(_sig, _kept)
        _w_m1      = rdistnn_weights_from_ranks(_rnk_m1, k, s, "asymp")
        _loo[_i,.] = cross(_w_m1, _sig_m1)
    }
    _centered = _loo :- mean(_loo)
    _V_bf = ((n - 1) / n) * cross(_centered, _centered)

    st_numscalar("_r_asymp_bf", max(abs(_V_fast - _V_bf)))
}
test_assert "[bf] asymp vector covariance matches brute-force LOO" scalar(_r_asymp_bf) < 1e-10

// ════════════════════════════════════════════════════════════
test_section "rdistnn_jk_cov_from_ranks: [shape] covariance properties"
// ════════════════════════════════════════════════════════════

mata {
    rseed(305)
    n = 18; q = 3; k = 2; s = 6
    _X   = rnormal(n, 2, 0, 1)
    _rnk = rdistnn_ranks(_X, (0.3, 0.1), "l2")
    _sig = rnormal(n, q, 0, 1)
    _V   = rdistnn_jk_cov_from_ranks(_rnk, _sig, k, s, "exact")

    _m1 = _V[1,1]
    _m2 = _V[2,2]
    _m3 = _V[3,3]
    _m12 = _V[1,1] * _V[2,2] - _V[1,2]^2
    _m13 = _V[1,1] * _V[3,3] - _V[1,3]^2
    _m23 = _V[2,2] * _V[3,3] - _V[2,3]^2
    _m123 = det(_V)
    _min_minor = min((_m1, _m2, _m3, _m12, _m13, _m23, _m123))

    st_numscalar("_r_sym", max(abs(_V - _V')))
    st_numscalar("_r_psd", _min_minor)
}
test_assert "[shape] vector covariance is symmetric" scalar(_r_sym) < 1e-10
test_assert "[shape] vector covariance is positive semidefinite" scalar(_r_psd) > -1e-8
