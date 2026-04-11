// _rdistnn_mata.ado
// Mata library for RdistNN: distance computation, rank ordering, DNN weight
// calculation, and rank-shift jackknife SE.  Loaded automatically when any
// rdistnn command is first called.
//
// Compile and save the .mlib via:
//   do dev-setup.do
//   run src/_rdistnn_mata.ado
//   quietly cd src
//   mata: mata mlib create lrdistnn, replace
//   mata: mata mlib add lrdistnn rdistnn_*()
//   quietly cd ..

// ---------------------------------------------------------------------
// Reload hygiene: if these functions already exist in the Mata session,
// drop them so this file can redefine the latest source version.
// ---------------------------------------------------------------------
capture mata: mata drop lncomb()
capture mata: mata drop rdistnn_l2dist()
capture mata: mata drop rdistnn_mahal_cov_inv()
capture mata: mata drop rdistnn_mahal_dist()
capture mata: mata drop rdistnn_order()
capture mata: mata drop rdistnn_weights_exact()
capture mata: mata drop rdistnn_weights_asymp()
capture mata: mata drop rdistnn_ranks()
capture mata: mata drop rdistnn_weights_from_ranks()
capture mata: mata drop rdistnn_weights()
capture mata: mata drop rdistnn_jk_se_from_ranks()

// ---------------------------------------------------------------------
// Plugin fast path: try to declare rdistnn_plugin once per session.
// Sets global $rdistnn_plugin_ok = 1 if the compiled plugin is found
// on adopath, 0 otherwise.  The main .ado commands check this global
// before falling back to the pure-Mata path.
// Note: the plugin only implements L2 distance.  Mahalanobis always
// uses the Mata path regardless of this global.
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

// ── lncomb ────────────────────────────────────────────────────────────
// Log binomial coefficient ln(C(n,k)) via lnfactorial().
// Returns -1e300 (effectively -inf) for out-of-support inputs.
real scalar lncomb(real scalar n, real scalar k)
{
    if (k < 0 | k > n | n < 0) return(-1e300)
    return(lnfactorial(n) - lnfactorial(k) - lnfactorial(n-k))
}

// ── rdistnn_l2dist ────────────────────────────────────────────────────
// Compute squared L2 distances from each row of X to the vector x0.
// Returns an n-vector of squared distances (unsorted).
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

// ── rdistnn_mahal_cov_inv ─────────────────────────────────────────────
// Compute the (generalized) inverse of the sample covariance matrix of X.
// Uses invsym() which applies the sweep operator and handles singular
// covariance gracefully (zeros out linearly dependent directions).
//
// Arguments:
//   X : n x kx matrix (kx >= 2)
//
// Returns: kx x kx inverse covariance matrix.
//
real matrix rdistnn_mahal_cov_inv(real matrix X)
{
    real matrix     Xc, S
    real scalar     n

    n  = rows(X)
    Xc = X :- mean(X)                // center each column
    S  = cross(Xc, Xc) / (n - 1)    // kx x kx sample covariance
    return(invsym(S))
}

// ── rdistnn_mahal_dist ────────────────────────────────────────────────
// Compute Mahalanobis distances from each row of X to the vector x0,
// given a precomputed inverse covariance matrix Sinv.
//
// d_i = sqrt( (x_i - x0)' * Sinv * (x_i - x0) )
//
// Arguments:
//   X    : n x kx matrix
//   x0   : 1 x kx evaluation point
//   Sinv : kx x kx inverse covariance (from rdistnn_mahal_cov_inv)
//
// Returns: n-vector of Mahalanobis distances (unsorted).
//
real colvector rdistnn_mahal_dist(real matrix X, real rowvector x0,
                                   real matrix Sinv)
{
    real matrix     Z
    real colvector  qform
    Z = X :- x0                        // broadcast row subtraction
    qform = rowsum((Z * Sinv) :* Z)
    // Clamp tiny negative roundoff to zero before sqrt().
    qform = qform :* (qform :>= 0)
    return(sqrt(qform))
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

// ── rdistnn_ranks ─────────────────────────────────────────────────────
// Compute the integer rank vector for X evaluated at x0, dispatching
// on the distance metric.  Encapsulates the L2 / Mahalanobis branch so
// that inline mata { } blocks in .ado files need no if-else braces.
//
// Arguments:
//   X        : n x kx matrix
//   x0       : 1 x kx evaluation point
//   distance : "l2" (default) or "mahalanobis"
//
// Returns: n-vector of integer ranks (1 = nearest).
//
real colvector rdistnn_ranks(real matrix X, real rowvector x0,
                              | string scalar distance)
{
    real matrix Sinv
    if (args() < 3 | distance == "" | distance == "l2") {
        return(rdistnn_order(rdistnn_l2dist(X, x0)))
    }
    else if (distance == "mahalanobis") {
        Sinv = rdistnn_mahal_cov_inv(X)
        return(rdistnn_order(rdistnn_mahal_dist(X, x0, Sinv)))
    }
    else {
        errprintf("rdistnn_ranks: unknown distance '%s'. Use 'l2' or 'mahalanobis'.\n", distance)
        exit(error(198))
    }
}

// ── rdistnn_weights_from_ranks ────────────────────────────────────────
// Compute DNN weights directly from a precomputed rank vector.
// Factored out of rdistnn_weights() so the jackknife can reuse it
// without recomputing distances.
//
// Arguments:
//   ranks  : n-vector of integer ranks (1 = nearest)
//   k      : number of nearest neighbors
//   s      : kernel order
//   method : "exact" (default) or "asymp"
//
// Returns: n-vector of DNN weights.
//
real colvector rdistnn_weights_from_ranks(real colvector ranks,
                                           real scalar k,
                                           real scalar s,
                                           | string scalar method)
{
    if (args() < 4 | method == "" | method == "exact") {
        return(rdistnn_weights_exact(ranks, k, s))
    }
    else if (method == "asymp") {
        return(rdistnn_weights_asymp(ranks, k, s))
    }
    else {
        errprintf("rdistnn_weights_from_ranks: unknown method '%s'. Use 'exact' or 'asymp'.\n", method)
        exit(error(198))
    }
}

// ── rdistnn_weights ───────────────────────────────────────────────────
// Top-level: compute DNN weights for X evaluated at x0.
//
// Arguments:
//   X        : n x kx matrix
//   x0       : 1 x kx evaluation point
//   k        : number of nearest neighbors
//   s        : kernel order
//   method   : "exact" (default) or "asymp"
//   distance : "l2" (default) or "mahalanobis"
//   Sinv     : precomputed inverse covariance (optional; computed if not supplied
//              and distance == "mahalanobis")
//
// Returns: n-vector of DNN weights.
//
real colvector rdistnn_weights(real matrix X, real rowvector x0,
                                real scalar k, real scalar s,
                                | string scalar method,
                                  string scalar distance,
                                  real matrix Sinv)
{
    real colvector dist, ranks

    if (args() < 5 | method == "")   method   = "exact"
    if (args() < 6 | distance == "") distance = "l2"

    if (distance == "mahalanobis") {
        if (args() < 7 | rows(Sinv) == 0) Sinv = rdistnn_mahal_cov_inv(X)
        dist = rdistnn_mahal_dist(X, x0, Sinv)
    }
    else {
        dist = rdistnn_l2dist(X, x0)
    }
    ranks = rdistnn_order(dist)
    return(rdistnn_weights_from_ranks(ranks, k, s, method))
}

// ── rdistnn_jk_se_from_ranks ──────────────────────────────────────────
// Rank-shift jackknife SE for a DNN-weighted sum.  When observation i
// (with rank r_i) is deleted, all remaining observations with rank > r_i
// are shifted down by 1, preserving the relative ordering under a fixed
// distance metric.
//
// This is exact delete-1 rank recomputation for L2.  For Mahalanobis, callers
// interpret the covariance metric as fixed at the full estimation sample.
// This matches the "fast" jackknife in the R reference (jackknife_rank_full).
//
// Arguments:
//   ranks   : n-vector of integer ranks from the full sample (1 = nearest)
//   signal  : n-vector of per-observation contributions (pseudo-outcomes, Y, ...)
//   k       : number of nearest neighbors
//   s       : kernel order
//   method  : "exact" (default) or "asymp"
//
// Returns: jackknife SE (scalar).
//
real scalar rdistnn_jk_se_from_ranks(real colvector ranks,
                                      real colvector signal,
                                      real scalar k,
                                      real scalar s,
                                      | string scalar method)
{
    real scalar     n, i, r_i
    real colvector  idx, kept, ranks_m1, signal_m1, w_m1, loo

    if (args() < 5 | method == "") method = "exact"

    n   = rows(ranks)
    loo = J(n, 1, .)
    idx = (1::n)

    for (i = 1; i <= n; i++) {
        r_i      = ranks[i]
        kept     = (idx :!= i)                            // n x 1 logical (0/1)
        ranks_m1 = select(ranks, kept)                    // n-1 ranks in 1..n with gap at r_i
        ranks_m1 = ranks_m1 :- (ranks_m1 :> r_i)         // shift down: now 1..n-1
        signal_m1 = select(signal, kept)
        w_m1     = rdistnn_weights_from_ranks(ranks_m1, k, s, method)
        loo[i]   = sum(w_m1 :* signal_m1)
    }

    return(sqrt((n-1) / n * sum((loo :- mean(loo)):^2)))
}

end
