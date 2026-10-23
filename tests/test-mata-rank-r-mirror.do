// test-mata-rank-r-mirror.do
// R-inspired deterministic rank, tie, and Mahalanobis edge tests.
//
// Mirrors the R package's ordering/distance batteries where the Stata Mata
// surface has direct equivalents.

quietly do "src/_rdistnn_mata.ado"
do "tests/test_framework.do"

test_section "rdistnn_order/rdistnn_ranks: [ties] deterministic tie handling"

mata {
    _d = (1 \ 1 \ 2 \ 2 \ 3 \ 3)
    _r = rdistnn_order(_d)
    st_numscalar("_tie_rank_diff", max(abs(_r - (1::6))))
}
test_assert "[ties] equal distances break by original observation order" scalar(_tie_rank_diff) == 0

mata {
    _X = (0,0 \ 0,0 \ 0,0 \ 0,0 \ 0,0)
    _r = rdistnn_ranks(_X, (0, 0), "l2")
    st_numscalar("_plateau_rank_diff", max(abs(_r - (1::5))))
}
test_assert "[ties] all-equal L2 plateau ranks as 1..n" scalar(_plateau_rank_diff) == 0

mata {
    _X = (1,0 \ -1,0 \ 0,1 \ 0,-1)
    _r1 = rdistnn_ranks(_X, (0, 0), "l2")
    _r2 = rdistnn_ranks(_X, (0, 0), "l2")
    st_numscalar("_tie_repeat", max(abs(_r1 - _r2)))
    st_numscalar("_tie_expected", max(abs(_r1 - (1::4))))
}
test_assert "[ties] repeated tie ranking is deterministic" scalar(_tie_repeat) == 0
test_assert "[ties] distinct symmetric L2 profiles use original row order" scalar(_tie_expected) == 0

test_section "rdistnn_ranks: [perm] row permutation maps ranks consistently"

mata {
    _X = (0,0 \ 2,0 \ 1,0 \ 4,0 \ 3,0)
    _p = (3 \ 1 \ 5 \ 2 \ 4)
    _r = rdistnn_ranks(_X, (0, 0), "l2")
    _rp = rdistnn_ranks(_X[_p,.], (0, 0), "l2")
    _mapped = J(rows(_X), 1, .)
    for (_i = 1; _i <= rows(_X); _i++) _mapped[_p[_i]] = _rp[_i]
    st_numscalar("_perm_rank_diff", max(abs(_r - _mapped)))
}
test_assert "[perm] ranks follow row permutation for unique distances" scalar(_perm_rank_diff) == 0

test_section "rdistnn_mahal_*: [edge] singular covariance and input constraints"

mata {
    _X = (1,3 \ 2,3 \ 3,3 \ 4,3 \ 5,3)
    _geometry = rdistnn_mahal_geometry(_X)
    _d = rdistnn_distance_matrix(_X, (3, 3), "mahalanobis", _geometry)
    _r = rdistnn_ranks(_X, (3, 3), "mahalanobis")
    st_numscalar("_mh_singular_finite", sum(_d :< .) == rows(_d))
    st_numscalar("_mh_singular_nonneg", min(_d) >= 0)
    st_numscalar("_mh_singular_perm", mreldif(sort(_r, 1), (1::rows(_r))))
}
test_assert "[edge] singular-cov Mahal distances are finite" scalar(_mh_singular_finite) == 1
test_assert "[edge] singular-cov Mahal distances are non-negative" scalar(_mh_singular_nonneg) == 1
test_assert "[edge] singular-cov Mahal ranks are a permutation" scalar(_mh_singular_perm) == 0

capture mata: _bad_ranks = rdistnn_ranks((1 \ 2 \ 3), (0), "cosine")
test_assert "[errors] unknown distance rejected by rdistnn_ranks()" _rc != 0
