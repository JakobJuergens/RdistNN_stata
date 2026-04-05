// _rdistnn_mata.ado
// Mata library for RdistNN: distance computation, rank ordering, and
// DNN weight calculation.  Loaded automatically when any rdistnn command
// is first called.
//
// Compile and save the .mlib via:
//   do dev-setup.do
//   run src/_rdistnn_mata.ado
//   quietly cd src
//   mata: mata mlib create lrdistnn, replace
//   mata: mata mlib add lrdistnn rdistnn_*()
//   quietly cd ..

// ---------------------------------------------------------------------
// Guard: only compile once per session
// ---------------------------------------------------------------------
capture mata: mata drop rdistnn_weights()
if _rc == 0 {
    // Already loaded; nothing to do.
    exit
}

mata:

// ── rdistnn_l2dist ────────────────────────────────────────────────────
// Compute squared L2 distances from each row of X to the vector x0.
// Returns an n-vector of distances (unsorted).
//
// Arguments:
//   X  : n x kx real matrix of covariate values
//   x0 : 1 x kx (or kx x 1) evaluation point
//
real colvector rdistnn_l2dist(real matrix X, real rowvector x0)
{
    real matrix  D
    D = X :- x0
    return(rowsum(D :* D))
}

// ── rdistnn_order ─────────────────────────────────────────────────────
// Return the integer rank order (1 = nearest) of each observation
// by distance to x0.  Ties broken by original index.
//
// Returns: n-vector of integer ranks.
//
real colvector rdistnn_order(real colvector dist2)
{
    real colvector  idx, rnk
    real scalar     n, i

    n   = rows(dist2)
    idx = order(dist2, 1)   // ascending order
    rnk = J(n, 1, .)
    for (i = 1; i <= n; i++) {
        rnk[idx[i]] = i
    }
    return(rnk)
}

// ── rdistnn_weights_exact ─────────────────────────────────────────────
// Compute exact DNN rank weights w_i(x0) = c(r_i, s, k) for each obs.
//
// The weight formula (from the paper, eq. for DNN rank weights):
//   w_i = sum_{j=r_i}^{s} C(j-1, k-1) / C(s, k)  for r_i <= s
//   w_i = 0                                          for r_i > s
//
// where r_i = rank of obs i by distance to x0,
//       k   = number of nearest neighbors,
//       s   = kernel order (bandwidth parameter).
//
// Arguments:
//   ranks : n-vector of integer ranks (1 = nearest)
//   k     : number of nearest neighbors (integer >= 1)
//   s     : kernel order (integer >= k)
//
// Returns: n-vector of weights summing to 1.
//
real colvector rdistnn_weights_exact(real colvector ranks,
                                      real scalar k,
                                      real scalar s)
{
    real scalar     n, i, r, j, num, denom
    real colvector  w

    n      = rows(ranks)
    w      = J(n, 1, 0)
    denom  = comb(s, k)

    if (denom == 0) {
        errprintf("rdistnn_weights_exact: s=%g < k=%g; invalid.\n", s, k)
        exit(error(198))
    }

    for (i = 1; i <= n; i++) {
        r = ranks[i]
        if (r > s) continue
        num = 0
        for (j = r; j <= s; j++) {
            num = num + comb(j - 1, k - 1)
        }
        w[i] = num / denom
    }
    return(w)
}

// ── rdistnn_weights ───────────────────────────────────────────────────
// Top-level: compute DNN weights for X evaluated at x0.
//
// Arguments:
//   X  : n x kx matrix
//   x0 : 1 x kx evaluation point
//   k  : number of nearest neighbors
//   s  : kernel order
//
// Returns: n-vector of DNN weights.
//
real colvector rdistnn_weights(real matrix X, real rowvector x0,
                                real scalar k, real scalar s)
{
    return(rdistnn_weights_exact(rdistnn_order(rdistnn_l2dist(X, x0)), k, s))
}

end
