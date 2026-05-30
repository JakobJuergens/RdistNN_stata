// test-mata-distances.do
// Tests for Mata distance and ranking functions:
//   rdistnn_l2dist, rdistnn_order, rdistnn_ranks,
//   rdistnn_mahal_cov_inv, rdistnn_mahal_dist
//
// Tags: [smoke] [edge] [inv] [known] [agree]
//
// Pattern: compute in a mata { } block, store result in a Stata scalar,
// then assert on the scalar with test_assert.

quietly do "src/_rdistnn_mata.ado"
do "tests/test_framework.do"

// ════════════════════════════════════════════════════════════
test_section "rdistnn_l2dist: [smoke] basic properties"
// ════════════════════════════════════════════════════════════

// [smoke] all L2 squared distances are non-negative
mata {
    rseed(1)
    _X = rnormal(50, 2, 0, 1)
    _d = rdistnn_l2dist(_X, (0, 0))
    st_numscalar("_r", min(_d))
}
test_assert "[smoke] L2 distances non-negative" scalar(_r) >= 0

// [edge] self-distance is exactly zero
mata {
    _X = (1, 2 \ 3, 4 \ 5, 6)
    _d = rdistnn_l2dist(_X, (3, 4))   // row 2 = evaluation point
    st_numscalar("_r", _d[2])
}
test_assert "[edge] L2 self-distance is 0" scalar(_r) == 0

// [edge] 1-D case: d[i] = (x[i] - x0)^2 exactly
mata {
    _X   = (1 \ 2 \ 3 \ 4 \ 5)
    _d   = rdistnn_l2dist(_X, (3))
    _exp = (4 \ 1 \ 0 \ 1 \ 4)
    st_numscalar("_r", mreldif(_d, _exp))
}
test_assert "[edge] 1D L2 distances match (x-x0)^2" scalar(_r) < 1e-12

// [inv] translation: adding same vector to X and x0 leaves distances unchanged
mata {
    rseed(2)
    _X  = rnormal(30, 2, 0, 1)
    _x0 = (0.5, -0.5)
    _d1 = rdistnn_l2dist(_X, _x0)
    _sh = (3.7, -2.1)
    _d2 = rdistnn_l2dist(_X :+ _sh, _x0 + _sh)
    st_numscalar("_r", mreldif(_d1, _d2))
}
test_assert "[inv] L2 translation-invariance" scalar(_r) < 1e-12

// ════════════════════════════════════════════════════════════
test_section "rdistnn_order: [smoke] permutation properties"
// ════════════════════════════════════════════════════════════

// [smoke] output is a permutation of 1..n
mata {
    rseed(3)
    _X   = rnormal(30, 2, 0, 1)
    _d   = rdistnn_l2dist(_X, (0, 0))
    _rnk = rdistnn_order(_d)
    _n   = rows(_rnk)
    st_numscalar("_r", mreldif(sort(_rnk, 1), (1::_n)))
}
test_assert "[smoke] rdistnn_order returns permutation of 1..n" scalar(_r) == 0

// [edge] nearest observation gets rank 1
mata {
    // obs at (3,0): nearest to x0=(3,0); 5 obs at x=1..5, y=0
    _X   = (1,0 \ 2,0 \ 3,0 \ 4,0 \ 5,0)
    _d   = rdistnn_l2dist(_X, (3, 0))
    _rnk = rdistnn_order(_d)
    st_numscalar("_r", _rnk[3])   // obs 3 should be rank 1
}
test_assert "[edge] nearest obs gets rank 1" scalar(_r) == 1

// [inv] rank order is invariant to monotone transformation of distances
//       (rdistnn_order on d vs sqrt(d) must agree)
mata {
    rseed(4)
    _X    = rnormal(40, 2, 0, 1)
    _x0   = (0.3, -0.7)
    _d2   = rdistnn_l2dist(_X, _x0)
    _d    = sqrt(_d2)
    _rnk2 = rdistnn_order(_d2)
    _rnkd = rdistnn_order(_d)
    st_numscalar("_r", max(abs(_rnk2 - _rnkd)))
}
test_assert "[inv] rank order invariant to sqrt of squared distances" scalar(_r) == 0

// ════════════════════════════════════════════════════════════
test_section "rdistnn_mahal_dist: [smoke] and [inv]"
// ════════════════════════════════════════════════════════════

// [smoke] Mahalanobis distances are finite and non-negative
mata {
    rseed(5)
    _X    = rnormal(50, 2, 0, 1)
    _x0   = (0, 0)
    _Sinv = rdistnn_mahal_cov_inv(_X)
    _d    = rdistnn_mahal_dist(_X, _x0, _Sinv)
    st_numscalar("_rmin", min(_d))
    st_numscalar("_rmax", max(_d))
    st_numscalar("_rfin", sum(_d :< .) == rows(_d))
}
test_assert "[smoke] Mahal distances non-negative" scalar(_rmin) >= 0
test_assert "[smoke] Mahal distances all finite"   scalar(_rfin) == 1

// [edge] Mahalanobis self-distance is (approximately) zero
mata {
    rseed(6)
    _X    = rnormal(40, 2, 0, 1)
    _Sinv = rdistnn_mahal_cov_inv(_X)
    // pick row 5 as evaluation point
    _d    = rdistnn_mahal_dist(_X, _X[5,], _Sinv)
    st_numscalar("_r", _d[5])
}
test_assert "[edge] Mahal self-distance is 0" scalar(_r) < 1e-10

// [inv] Mahalanobis distance is scale-invariant:
//       if Y = c*X and y0 = c*x0, then Mahal_Y(yi,y0) = Mahal_X(xi,x0)
mata {
    rseed(7)
    _X    = rnormal(50, 2, 0, 1)
    _x0   = (0.4, -0.3)
    _Sinv1 = rdistnn_mahal_cov_inv(_X)
    _d1    = rdistnn_mahal_dist(_X, _x0, _Sinv1)
    // scale by 3
    _Y     = 3 * _X
    _y0    = 3 * _x0
    _Sinv2 = rdistnn_mahal_cov_inv(_Y)
    _d2    = rdistnn_mahal_dist(_Y, _y0, _Sinv2)
    st_numscalar("_r", mreldif(_d1, _d2))
}
test_assert "[inv] Mahal distances scale-invariant" scalar(_r) < 1e-10

// [known] diagonal covariance: Mahal_i = sqrt(x1i^2/s1 + x2i^2/s2)
//   Construct X so that sample cov = diag(4, 1) exactly (n=5, Bessel-corrected).
//   col1 = (2, 2,-2,-2, 0)': SS=16, var=4  ✓
//   col2 = (1,-1, 1,-1, 0)': SS= 4, var=1  ✓
mata {
    _X  = (2,1 \ 2,-1 \ -2,1 \ -2,-1 \ 0,0)
    _x0 = (0, 0)
    _Sinv = rdistnn_mahal_cov_inv(_X)
    _d    = rdistnn_mahal_dist(_X, _x0, _Sinv)
    // expected: sqrt(x1^2/4 + x2^2/1)
    _exp  = sqrt(_X[,1]:^2 / 4 + _X[,2]:^2 / 1)
    st_numscalar("_r", mreldif(_d, _exp))
}
test_assert "[known] diagonal-cov Mahal matches closed form" scalar(_r) < 1e-10

// ════════════════════════════════════════════════════════════
test_section "rdistnn_ranks: [agree] dispatch"
// ════════════════════════════════════════════════════════════

// [agree] rdistnn_ranks("l2") matches rdistnn_order(rdistnn_l2dist(...))
mata {
    rseed(8)
    _X   = rnormal(40, 2, 0, 1)
    _x0  = (0.1, -0.2)
    _r1  = rdistnn_ranks(_X, _x0, "l2")
    _r2  = rdistnn_order(rdistnn_l2dist(_X, _x0))
    st_numscalar("_r", max(abs(_r1 - _r2)))
}
test_assert "[agree] rdistnn_ranks l2 = order(l2dist)" scalar(_r) == 0

// [agree] rdistnn_ranks("mahalanobis") matches manual Mahal pipeline
mata {
    rseed(9)
    _X    = rnormal(40, 2, 0, 1)
    _x0   = (0.2, 0.1)
    _r1   = rdistnn_ranks(_X, _x0, "mahalanobis")
    _Sinv = rdistnn_mahal_cov_inv(_X)
    _r2   = rdistnn_order(rdistnn_mahal_dist(_X, _x0, _Sinv))
    st_numscalar("_r", max(abs(_r1 - _r2)))
}
test_assert "[agree] rdistnn_ranks mahal = order(mahal_dist)" scalar(_r) == 0
