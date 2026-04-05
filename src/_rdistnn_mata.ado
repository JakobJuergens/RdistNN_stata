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

// ---------------------------------------------------------------------
// Plugin fast path: try to declare rdistnn_plugin once per session.
// Sets global $rdistnn_plugin_ok = 1 if the compiled plugin is found
// on adopath, 0 otherwise.  The main .ado commands check this global
// before falling back to the pure-Mata path.
// ---------------------------------------------------------------------
if ("$rdistnn_plugin_ok" == "") {
    // _rc == 0   : declared for the first time
    // _rc == 110 : already declared in this session
    // anything else: .plugin file not found on adopath
    capture program define rdistnn_plugin, plugin
    if _rc == 0 | _rc == 110 {
        global rdistnn_plugin_ok 1
    }
    else {
        global rdistnn_plugin_ok 0
    }
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
// Exact DNN rank weights (Steele closed form / full WOR hypergeometric).
//
// For obs with global rank r (1 = nearest to x0):
//
//   w_r = (1/k) * sum_{l=1}^{k} C(r-1, l-1) * C(n-r, s-l) / C(n, s)
//
// Support: r = 1, ..., n-s+k  (NOT 1..s).  For k=1 this equals
//   C(n-r, s-1) / C(n, s), which is the standard Steele formula.
//
// Derivation: draw a uniformly random s-subset from all n observations;
// find the k nearest to x0 within that subset.  The weight for obs i is
// the marginal probability that it is among those k nearest, summed over
// all C(n,s) equally-likely subsets.
//
// Arguments:
//   ranks : n-vector of integer ranks (1 = nearest); n = rows(ranks)
//   k     : number of nearest neighbors (integer >= 1)
//   s     : kernel order / subset size (integer, k <= s <= n)
//
// Returns: n-vector of weights summing to 1 (up to floating-point error).
//
real colvector rdistnn_weights_exact(real colvector ranks,
                                      real scalar k,
                                      real scalar s)
{
    real scalar     n, i, r, l, l_lo, l_hi, lp, log_C_n_s, wsum
    real colvector  w

    n = rows(ranks)
    w = J(n, 1, 0)

    if (s < k | s > n) {
        errprintf("rdistnn_weights_exact: need k <= s <= n; got k=%g s=%g n=%g\n", k, s, n)
        exit(error(198))
    }

    log_C_n_s = lncomb(n, s)

    for (i = 1; i <= n; i++) {
        r = ranks[i]
        // Support: r <= n - s + k
        if (r > n - s + k) continue

        // l must satisfy: l >= 1, l <= min(k, r, s),
        //                 and l >= s - (n - r)  [so C(n-r, s-l) is nonneg]
        l_lo = max((1, s - (n - r)))
        l_hi = min((k, min((r, s))))
        if (l_lo > l_hi) continue

        lp = 0
        for (l = l_lo; l <= l_hi; l++) {
            lp = lp + exp(lncomb(r-1, l-1) + lncomb(n-r, s-l) - log_C_n_s)
        }
        w[i] = lp / k
    }

    // Normalize for floating-point safety (should sum to 1 analytically)
    wsum = sum(w)
    if (wsum > 0) w = w / wsum
    return(w)
}

// ── rdistnn_weights_asymp ─────────────────────────────────────────────
// Asymptotic DNN weights via Beta order-statistic increments.
//
// For rank j (1-based), the asymptotic weight is:
//   w_rank[j] = (1/k) * sum_{l=1}^{k} [ F_Beta(l, s-l+1)(j/n)
//                                       - F_Beta(l, s-l+1)((j-1)/n) ]
//
// This is the large-n approximation to the exact WOR formula, treating
// the sorted s-subset positions as i.i.d. Uniform(0,1) order statistics.
// Full support over 1..n; renormalized to sum exactly 1.
//
// Arguments:
//   ranks : n-vector of integer ranks (1 = nearest)
//   k     : number of nearest neighbors (integer >= 1)
//   s     : kernel order (integer >= k)
//
// Returns: n-vector of weights summing to 1.
//
real colvector rdistnn_weights_asymp(real colvector ranks,
                                      real scalar k,
                                      real scalar s)
{
    real scalar     n, i, l, r
    real colvector  w, w_rank

    n      = rows(ranks)
    w_rank = J(n, 1, 0)

    for (l = 1; l <= k; l++) {
        // Beta(l, s-l+1) CDF increments on ((j-1)/n, j/n] for j=1..n
        for (i = 1; i <= n; i++) {
            w_rank[i] = w_rank[i] + (ibeta(l, s-l+1, i/n) - ibeta(l, s-l+1, (i-1)/n))
        }
    }
    w_rank = w_rank / k

    // Normalize for safety
    w_rank = w_rank / sum(w_rank)

    // Map from rank order back to observation order
    w = J(n, 1, 0)
    for (i = 1; i <= n; i++) {
        r    = ranks[i]
        w[i] = w_rank[r]
    }
    return(w)
}

// ── rdistnn_weights ───────────────────────────────────────────────────
// Top-level: compute DNN weights for X evaluated at x0.
//
// Arguments:
//   X      : n x kx matrix
//   x0     : 1 x kx evaluation point
//   k      : number of nearest neighbors
//   s      : kernel order
//   method : "exact" (default) or "asymp"
//
// Returns: n-vector of DNN weights.
//
real colvector rdistnn_weights(real matrix X, real rowvector x0,
                                real scalar k, real scalar s,
                                | string scalar method)
{
    real colvector ranks
    if (args() < 5 | method == "" | method == "exact") {
        ranks = rdistnn_order(rdistnn_l2dist(X, x0))
        return(rdistnn_weights_exact(ranks, k, s))
    }
    else if (method == "asymp") {
        ranks = rdistnn_order(rdistnn_l2dist(X, x0))
        return(rdistnn_weights_asymp(ranks, k, s))
    }
    else {
        errprintf("rdistnn_weights: unknown method '%s'. Use 'exact' or 'asymp'.\n", method)
        exit(error(198))
    }
}

end
