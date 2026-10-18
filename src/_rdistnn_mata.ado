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
capture mata: mata drop rdistnn_distance_matrix()
capture mata: mata drop rdistnn_ranks_from_dmat()
capture mata: mata drop rdistnn_ranks_matrix()
capture mata: mata drop rdistnn_weighted_quantile()
capture mata: mata drop rdistnn_weights_from_ranks()
capture mata: mata drop rdistnn_weights_from_rmat()
capture mata: mata drop rdistnn_weights()
capture mata: mata drop rdistnn_loo_naive()
capture mata: mata drop rdistnn_loo_prefix()
capture mata: mata drop rdistnn_loo_plugin()
capture mata: mata drop rdistnn_loo_from_ranks()
capture mata: mata drop rdistnn_center_signal()
capture mata: mata drop rdistnn_jk_cov_from_ranks()
capture mata: mata drop rdistnn_jk_cov_multi_from_ranks()
capture mata: mata drop rdistnn_jk_se_from_ranks()
capture mata: mata drop rdistnn_plot_target_diag()
capture mata: mata drop rdistnn_plot_candidate_diag()
capture mata: mata drop rdistnn_plot_est_default()
capture mata: mata drop rdistnn_plot_est_grid()
capture mata: mata drop rdistnn_plot_tuning()
capture mata: mata drop rdistnn_plot_overlap()
capture mata: mata drop rdistnn_plot_rank_weights()
capture mata: mata drop rdistnn_plot_kernel()
capture mata: mata drop rdistnn_plot_influence()

// ---------------------------------------------------------------------
// Plugin fast path: try to declare rdistnn_plugin once per session.
// Sets global $rdistnn_plugin_ok = 1 if the compiled plugin is found
// on adopath, 0 otherwise.  Commands may still choose the pure-Mata path
// unless plugin use is explicitly enabled for parity testing.
// Note: the plugin only implements L2 distance.  Mahalanobis always
// uses the Mata path regardless of this global.
// ---------------------------------------------------------------------
if ("$rdistnn_plugin_ok" == "") {
    // _rc == 0   : declared for the first time
    // _rc == 110 : already declared in this session
    // anything else: rdistnn.plugin file not found on adopath
    capture program define rdistnn_plugin, plugin using(rdistnn.plugin)
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
// Match R's ordinary inverse with an SVD Moore-Penrose fallback. A sweep
// generalized inverse can change distances for targets outside the data span.
//
// Arguments:
//   X : n x kx matrix (kx >= 2)
//
// Returns: kx x kx inverse covariance matrix.
//
real matrix rdistnn_mahal_cov_inv(real matrix X)
{
    real matrix     Xc, S, Si, U, Vt
    real colvector  values, inverse_values
    real scalar     n, condition_product, cutoff, j

    n  = rows(X)
    if (n < 2) {
        errprintf("rdistnn_mahal_cov_inv: Mahalanobis distance requires at least two complete observations.\n")
        exit(error(2001))
    }
    Xc = X :- mean(X)                // center each column
    S  = cross(Xc, Xc) / (n - 1)    // kx x kx sample covariance
    Si = luinv(S, 0)
    if (!any(Si :>= .)) {
        condition_product = norm(S, 1) * norm(Si, 1)
        // R's solve() rejects reciprocal condition numbers below double eps.
        if (condition_product < . & condition_product > 0) {
            if (1/condition_product >= 2.220446049250313e-16) return(Si)
        }
    }
    svd(S, U, values, Vt)
    cutoff = 1e-10 * max((max(values), 1e-300))
    inverse_values = J(rows(values), 1, 0)
    for (j=1; j<=rows(values); j++) {
        if (values[j] > cutoff) inverse_values[j] = 1/values[j]
    }
    return(Vt' * diag(inverse_values) * U')
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
    idx = order((dist2, (1::n)), (1, 2))   // ascending distance, then row index
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
// Positive-weight support: r = 1, ..., n-s+k  (NOT 1..s).  Thus larger s
// shrinks support and tightens localization; for k=1 the edge is n-s+1.
// Here k is the number selected within each s-subset, not the global support
// size.  For k=1 the weight equals C(n-r, s-1) / C(n, s), the Steele formula.
//
// Arguments:
//   ranks : n-vector of integer ranks (1 = nearest); n = rows(ranks)
//   k     : number selected within each s-subset (integer >= 1)
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
//   k     : number selected within each s-subset (integer >= 1)
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

// Compute an n x m matrix of ordinary distances for all target rows. L2
// ranking may use squared distances, but user-facing diagnostics use the
// Euclidean scale and therefore take the square root here.
real matrix rdistnn_distance_matrix(real matrix X, real matrix X0,
                                     | string scalar distance)
{
    real scalar j, m
    real matrix dmat, Sinv

    if (cols(X0) != cols(X)) {
        errprintf("rdistnn_distance_matrix: cols(X0) must equal cols(X).\n")
        exit(error(198))
    }
    m = rows(X0)
    if (m < 1) {
        errprintf("rdistnn_distance_matrix: X0 must contain at least one target row.\n")
        exit(error(198))
    }

    dmat = J(rows(X), m, .)
    if (args() < 3 | distance == "" | distance == "l2") {
        for (j = 1; j <= m; j++) {
            dmat[,j] = sqrt(rdistnn_l2dist(X, X0[j,.]))
        }
    }
    else if (distance == "mahalanobis") {
        Sinv = rdistnn_mahal_cov_inv(X)
        for (j = 1; j <= m; j++) {
            dmat[,j] = rdistnn_mahal_dist(X, X0[j,.], Sinv)
        }
    }
    else {
        errprintf("rdistnn_distance_matrix: unknown distance '%s'. Use 'l2' or 'mahalanobis'.\n", distance)
        exit(error(198))
    }
    return(dmat)
}

// Convert each distance column to a full rank permutation.
real matrix rdistnn_ranks_from_dmat(real matrix dmat)
{
    real scalar j
    real matrix ranks

    if (rows(dmat) < 1 | cols(dmat) < 1) {
        errprintf("rdistnn_ranks_from_dmat: dmat must be nonempty.\n")
        exit(error(198))
    }
    ranks = J(rows(dmat), cols(dmat), .)
    for (j = 1; j <= cols(dmat); j++) {
        ranks[,j] = rdistnn_order(dmat[,j])
    }
    return(ranks)
}

// Rank every observation separately for each row of X0. The returned n x m
// matrix has one full rank permutation per target column.
real matrix rdistnn_ranks_matrix(real matrix X, real matrix X0,
                                  | string scalar distance)
{
    if (args() < 3) distance = "l2"
    return(rdistnn_ranks_from_dmat(
        rdistnn_distance_matrix(X, X0, distance)))
}

// Left-continuous weighted empirical quantile used by R diagnostics.
real scalar rdistnn_weighted_quantile(real colvector x,
                                       real colvector w,
                                       real scalar prob)
{
    real scalar i, accum, compensation, adjusted, next_accum, wsum
    real colvector ord

    if (rows(x) != rows(w) | rows(x) < 1) {
        errprintf("rdistnn_weighted_quantile: x and w must have the same positive length.\n")
        exit(error(198))
    }
    if (any(missing(x))) {
        errprintf("rdistnn_weighted_quantile: x must contain only finite, nonmissing values.\n")
        exit(error(198))
    }
    if (missing(prob) | prob < 0 | prob > 1) {
        errprintf("rdistnn_weighted_quantile: probability must be finite and lie in [0,1].\n")
        exit(error(198))
    }
    if (any(missing(w)) | any(w :< 0)) {
        errprintf("rdistnn_weighted_quantile: weights must be finite and nonnegative.\n")
        exit(error(198))
    }
    wsum = quadsum(w)
    if (missing(wsum) | wsum <= 0) {
        errprintf("rdistnn_weighted_quantile: total weight must be strictly positive.\n")
        exit(error(198))
    }
    ord = order((x, (1::rows(x))), (1, 2))
    accum = 0
    compensation = 0
    for (i = 1; i <= rows(ord); i++) {
        adjusted = w[ord[i]] - compensation
        next_accum = accum + adjusted
        compensation = (next_accum - accum) - adjusted
        accum = next_accum
        if (accum / wsum >= prob) return(x[ord[i]])
    }
    return(x[ord[rows(ord)]])
}

// ── rdistnn_weights_from_ranks ────────────────────────────────────────
// Compute DNN weights directly from a precomputed rank vector.
// Factored out of rdistnn_weights() so the jackknife can reuse it
// without recomputing distances.
//
// Arguments:
//   ranks  : n-vector of integer ranks (1 = nearest)
//   k      : number selected within each s-subset
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

// Apply the rank-to-weight map independently to every target/rank column.
real matrix rdistnn_weights_from_rmat(real matrix ranks,
                                       real scalar k,
                                       real scalar s,
                                       | string scalar method)
{
    real scalar j, m
    real colvector table
    real matrix W

    if (args() < 4 | method == "") method = "exact"
    m = cols(ranks)
    if (m < 1) {
        errprintf("rdistnn_weights_from_rmat: ranks must have at least one column.\n")
        exit(error(198))
    }

    // The rank-weight table is independent of the target permutation.
    W = J(rows(ranks), m, .)
    table = rdistnn_weights_from_ranks((1::rows(ranks)), k, s, method)
    for (j = 1; j <= m; j++) {
        W[,j] = table[ranks[,j]]
    }
    return(W)
}

// ── rdistnn_weights ───────────────────────────────────────────────────
// Top-level: compute DNN weights for X evaluated at x0.
//
// Arguments:
//   X        : n x kx matrix
//   x0       : 1 x kx evaluation point
//   k        : number selected within each s-subset
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

// Explicit rank-shift deletions retained as a reference and final fallback.
// Output columns are target blocks, with q signal columns in each block.
real matrix rdistnn_loo_naive(real matrix ranks, real matrix signal,
    real scalar k, real scalar s, string scalar method)
{
    real scalar n, m, q, i, j
    real colvector kept, r, w
    real matrix loo
    n=rows(ranks); m=cols(ranks); q=cols(signal)
    loo=J(n,m*q,.)
    for(i=1;i<=n;i++) {
        kept=((1::n):!=i)
        for(j=1;j<=m;j++) {
            r=select(ranks[,j],kept)
            r=r:-(r:>ranks[i,j])
            w=rdistnn_weights_from_ranks(r,k,s,method)
            loo[i,((j-1)*q+1)..(j*q)]=cross(w,select(signal,kept))
        }
    }
    return(loo)
}

// Two linear passes per target/component. The table at n-1 is shared across
// every deletion and target; no deletion-specific rank or signal copies.
real matrix rdistnn_loo_prefix(real matrix ranks, real matrix signal,
    real colvector table)
{
    real scalar n, m, q, j, r
    real colvector ord
    real matrix z, left, right, loo
    n=rows(ranks); m=cols(ranks); q=cols(signal)
    loo=J(n,m*q,.)
    for(j=1;j<=m;j++) {
        ord=J(n,1,.)
        ord[ranks[,j]]=(1::n)
        z=signal[ord,]
        left=right=J(n,q,0)
        for(r=2;r<=n;r++) left[r,]=left[r-1,]+table[r-1]*z[r-1,]
        for(r=n-1;r>=1;r--) right[r,]=right[r+1,]+table[r]*z[r+1,]
        loo[,((j-1)*q+1)..(j*q)]=left[ranks[,j],]+right[ranks[,j],]
    }
    return(loo)
}

// The bridge uses Mata-created Stata matrices, which may exceed c(matsize).
// Attempt the native engine without changing the user's matrix setting.
real matrix rdistnn_loo_plugin(real matrix ranks, real matrix signal,
    real colvector table)
{
    string rowvector names
    real scalar rc, i
    real matrix loo
    if(st_global("rdistnn_plugin_ok")!="1") return(J(0,0,.))
    names=st_tempname(4)
    st_matrix(names[1],ranks)
    st_matrix(names[2],signal)
    st_matrix(names[3],table)
    st_matrix(names[4],J(rows(ranks),cols(ranks)*cols(signal),.))
    rc=_stata("plugin call rdistnn_plugin, rank_loo "+invtokens(names),1)
    loo=J(0,0,.)
    if(rc==0) loo=st_matrix(names[4])
    for(i=1;i<=4;i++) st_matrix(names[i],J(0,0,.))
    return(loo)
}

// Default hierarchy: compiled prefix/suffix -> Mata prefix/suffix -> naive.
// Explicit engine selectors support parity checks without changing estimators.
real matrix rdistnn_loo_from_ranks(real matrix ranks, real matrix signal,
    real scalar k, real scalar s, | string scalar method)
{
    real scalar n, m, q, i, j, r
    real colvector seen, table
    real matrix loo
    string scalar engine, plugin_mode
    pointer(function) scalar prefix
    if(args()<5 | method=="") method="exact"
    n=rows(ranks); m=cols(ranks); q=cols(signal)
    if((method!="exact" & method!="asymp") |
        n<2 | m<1 | rows(signal)!=n | q<1 | hasmissing(signal) |
        hasmissing(ranks) | k<1 | k!=floor(k) | s<k | s>=n | s!=floor(s) |
        min(ranks)<1 | max(ranks)>n | any(ranks:!=floor(ranks))) {
        errprintf("rdistnn_loo_from_ranks: require finite aligned signals, rank permutations and 1 <= k <= s < n.\n")
        exit(error(198))
    }
    for(j=1;j<=m;j++) {
        seen=J(n,1,0)
        for(i=1;i<=n;i++) {
            r=ranks[i,j]
            if(seen[r]) {
                errprintf("rdistnn_loo_from_ranks: each ranks column must be a permutation of 1..n.\n")
                exit(error(198))
            }
            seen[r]=1
        }
    }
    engine=st_global("rdistnn_jackknife_engine")
    if(engine=="") engine="auto"
    if(!any(engine:==("auto","plugin","mata","naive"))) {
        errprintf("rdistnn_jackknife_engine must be auto, plugin, mata or naive.\n")
        exit(error(198))
    }
    if(engine=="naive") return(rdistnn_loo_naive(ranks,signal,k,s,method))
    // CLATE delegates geometry/tuning with a zero signal. Its unused scalar
    // covariance needs neither a weight table nor any deletion calculations.
    if(all(signal:==0)) return(J(n,m*q,0))
    table=rdistnn_weights_from_ranks((1::(n-1)),k,s,method)
    plugin_mode=st_global("rdistnn_use_plugin")
    if(engine=="plugin" | (engine=="auto" & plugin_mode!="0")) {
        loo=rdistnn_loo_plugin(ranks,signal,table)
        if(rows(loo)==n & !hasmissing(loo)) return(loo)
        if(engine=="plugin" | plugin_mode=="1") {
            errprintf("Compiled rank_loo unavailable or failed; install an updated plugin or allow Mata fallback.\n")
            exit(error(498))
        }
    }
    // Resolve dynamically so the retained reference still works if the fast
    // Mata function is unavailable. Nonfinite fast sums also retry the reference.
    prefix=findexternal("rdistnn_loo_prefix()")
    if(prefix!=NULL) {
        loo=(*prefix)(ranks,signal,table)
        if(!hasmissing(loo)) return(loo)
    }
    if(engine=="mata") {
        errprintf("Mata prefix/suffix engine unavailable or produced nonfinite deletions.\n")
        exit(error(498))
    }
    return(rdistnn_loo_naive(ranks,signal,k,s,method))
}

// A common signal level has zero jackknife covariance. Remove it before
// accumulating deletions, and never restore it on covariance-only paths.
// Use the same-sign endpoint nearest zero; mixed-sign columns need no shift.
// A midpoint or a remote outlier would erase small, locally relevant signals.
real matrix rdistnn_center_signal(real matrix signal)
{
    real scalar j
    real rowvector lo, hi, level
    lo=colmin(signal); hi=colmax(signal)
    level=J(1,cols(signal),0)
    for(j=1;j<=cols(signal);j++) {
        if(lo[j]>0) level[j]=lo[j]
        else if(hi[j]<0) level[j]=hi[j]
    }
    return(signal:-level)
}

// ── rdistnn_jk_cov_from_ranks ─────────────────────────────────────────
// Rank-shift jackknife covariance for a DNN-weighted signal matrix.  When
// observation i
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
//   signal  : n x q matrix of per-observation contributions
//   k       : number selected within each s-subset
//   s       : kernel order
//   method  : "exact" (default) or "asymp"
//
// Returns: q x q jackknife covariance matrix.
//
real matrix rdistnn_jk_cov_from_ranks(real colvector ranks,
                                       real matrix signal,
                                       real scalar k,
                                       real scalar s,
                                       | string scalar method)
{
    real scalar     n, q
    real matrix     loo, centered, V

    if (args() < 5 | method == "") method = "exact"

    n = rows(ranks)
    if (rows(signal) != n) {
        errprintf("rdistnn_jk_cov_from_ranks: rows(signal) must equal rows(ranks).\n")
        exit(error(198))
    }

    q = cols(signal)
    if (q < 1) {
        errprintf("rdistnn_jk_cov_from_ranks: signal must have at least one column.\n")
        exit(error(198))
    }

    loo = rdistnn_loo_from_ranks(ranks, rdistnn_center_signal(signal), k, s, method)
    if(all(loo:==0)) return(J(q,q,0))

    centered = loo :- mean(loo)
    V = ((n - 1) / n) * cross(centered, centered)
    return((V + V') / 2)
}

// Joint rank-shift jackknife covariance for a scalar signal evaluated at
// multiple targets. Each ranks column is the full-sample rank permutation for
// one target. The returned m x m matrix retains cross-target covariance.
real matrix rdistnn_jk_cov_multi_from_ranks(real matrix ranks,
                                             real colvector signal,
                                             real scalar k,
                                             real scalar s,
                                             | string scalar method)
{
    real scalar     n, m
    real matrix     loo, centered, V

    if (args() < 5 | method == "") method = "exact"
    n = rows(ranks)
    m = cols(ranks)
    if (rows(signal) != n | cols(signal) != 1) {
        errprintf("rdistnn_jk_cov_multi_from_ranks: signal must be an n x 1 vector matching ranks.\n")
        exit(error(198))
    }
    if (m < 1) {
        errprintf("rdistnn_jk_cov_multi_from_ranks: ranks must have at least one target column.\n")
        exit(error(198))
    }

    loo = rdistnn_loo_from_ranks(ranks, rdistnn_center_signal(signal), k, s, method)
    if(all(loo:==0)) return(J(m,m,0))

    centered = loo :- mean(loo)
    V = ((n - 1) / n) * cross(centered, centered)
    return((V + V') / 2)
}

// ── rdistnn_jk_se_from_ranks ──────────────────────────────────────────
// Scalar compatibility wrapper for callers that only need the marginal
// jackknife SE of one signal column.
//
real scalar rdistnn_jk_se_from_ranks(real colvector ranks,
                                      real colvector signal,
                                      real scalar k,
                                      real scalar s,
                                      | string scalar method)
{
    real matrix V

    if (args() < 5 | method == "") method = "exact"
    V = rdistnn_jk_cov_from_ranks(ranks, signal, k, s, method)
    return(sqrt(max((V[1,1], 0))))
}

// ── postestimation plot-data builders ────────────────────────────────
// These helpers keep Mata loops outside ado if/else braces, whose parser can
// otherwise detach the matching else branch. They return plain numeric
// matrices consumed by cate_dnn_diagnostics and cate_dnn_plot.

real matrix rdistnn_plot_target_diag(real matrix diagnostics,
                                              real matrix targets)
{
    if (rows(diagnostics) != rows(targets)) {
        errprintf("rdistnn_plot_target_diagnostics: row counts must match.\n")
        exit(error(503))
    }
    return(((1::rows(diagnostics)), targets, diagnostics))
}

real matrix rdistnn_plot_candidate_diag(real matrix meta,
                                                 real matrix diagnostics,
                                                 real colvector idx)
{
    if (rows(meta) != rows(diagnostics)) {
        errprintf("rdistnn_plot_candidate_diagnostics: row counts must match.\n")
        exit(error(503))
    }
    return((meta[idx,], diagnostics[idx,]))
}

real matrix rdistnn_plot_est_default(real matrix targets,
                                            real colvector estimate,
                                            real colvector se,
                                            real scalar s,
                                            real scalar z)
{
    real scalar m
    m = rows(targets)
    if (rows(estimate) != m | rows(se) != m) {
        errprintf("rdistnn_plot_estimates_default: estimate and se must match targets.\n")
        exit(error(503))
    }
    return((targets, J(m, 1, 0), J(m, 1, s), J(m, 1, 1),
        (1::m), estimate, se, estimate :- z :* se, estimate :+ z :* se))
}

real matrix rdistnn_plot_est_grid(real matrix targets,
                                         real matrix meta,
                                         real matrix estimates,
                                         real matrix ses,
                                         real colvector idx,
                                         real scalar z)
{
    real scalar i, g, m
    real colvector b, se
    real matrix out

    m = rows(targets)
    out = J(0, cols(targets) + 8, .)
    for (i = 1; i <= rows(idx); i++) {
        g = idx[i]
        b = estimates[g,]'
        se = ses[g,]'
        out = out \ (targets, J(m, 1, g), J(m, 1, meta[g,8]),
            J(m, 1, meta[g,3]), (1::m), b, se,
            b :- z :* se, b :+ z :* se)
    }
    return(out)
}

real matrix rdistnn_plot_tuning(real matrix meta,
                                      real matrix estimates,
                                      real matrix ses,
                                      real colvector candidate_idx,
                                      real colvector target_idx,
                                      real scalar z)
{
    real scalar i, j, g, t, b, se
    real matrix out
    out = J(0, 9, .)
    for (i = 1; i <= rows(candidate_idx); i++) {
        g = candidate_idx[i]
        for (j = 1; j <= rows(target_idx); j++) {
            t = target_idx[j]
            b = estimates[g,t]
            se = ses[g,t]
            out = out \ (g, meta[g,8], meta[g,7], meta[g,3], t,
                b, se, b - z * se, b + z * se)
        }
    }
    return(out)
}

real matrix rdistnn_plot_overlap(real matrix diagnostics,
                                       real matrix targets,
                                       real rowvector cols_use)
{
    if (rows(diagnostics) != rows(targets)) {
        errprintf("rdistnn_plot_overlap_data: row counts must match.\n")
        exit(error(503))
    }
    return((targets, (1::rows(diagnostics)), diagnostics[,cols_use]))
}

real matrix rdistnn_plot_rank_weights(real matrix meta,
                                       real colvector idx,
                                       real scalar n,
                                       real scalar k,
                                       real scalar default_s,
                                       string scalar method,
                                       real scalar use_meta)
{
    real scalar i, g, s, is_default
    real colvector w
    real matrix out
    out = J(0, 6, .)
    for (i = 1; i <= rows(idx); i++) {
        g = idx[i]
        if (use_meta) {
            s = meta[g,8]
            is_default = meta[g,3]
        }
        else {
            s = default_s
            is_default = 1
        }
        w = rdistnn_weights_from_ranks((1::n), k, s, method)
        out = out \ ((1::n), w, runningsum(w), J(n,1,g),
            J(n,1,s), J(n,1,is_default))
    }
    return(out)
}

real matrix rdistnn_plot_kernel(real matrix X,
                                      real matrix targets,
                                      real colvector target_idx,
                                      real scalar k,
                                      real scalar s,
                                      string scalar method,
                                      string scalar distance)
{
    real scalar j, t, n
    real colvector ord, cum
    real matrix DX, ranks, W, out
    n = rows(X)
    DX = rdistnn_distance_matrix(X, targets, distance)
    ranks = rdistnn_ranks_from_dmat(DX)
    W = rdistnn_weights_from_rmat(ranks, k, s, method)
    out = J(0, 5, .)
    for (j = 1; j <= rows(target_idx); j++) {
        t = target_idx[j]
        ord = order((ranks[,t], (1::n)), (1,2))
        cum = J(n,1,.)
        cum[ord] = runningsum(W[ord,t])
        out = out \ (J(n,1,t), ranks[,t], DX[,t], W[,t], cum)
    }
    return(out)
}

real matrix rdistnn_plot_influence(real matrix X,
                                         real matrix targets,
                                         real colvector rho,
                                         real colvector obs_id,
                                         real rowvector estimates,
                                         real colvector target_idx,
                                         real scalar k,
                                         real scalar s,
                                         string scalar method,
                                         string scalar distance)
{
    real scalar j, t
    real colvector contribution
    real matrix DX, ranks, W, out
    DX = rdistnn_distance_matrix(X, targets, distance)
    ranks = rdistnn_ranks_from_dmat(DX)
    W = rdistnn_weights_from_rmat(ranks, k, s, method)
    out = J(0, 7, .)
    for (j = 1; j <= rows(target_idx); j++) {
        t = target_idx[j]
        contribution = W[,t] :* (rho :- estimates[1,t])
        out = out \ (J(rows(X),1,t), obs_id, ranks[,t], DX[,t],
            W[,t], rho, contribution)
    }
    return(out)
}

end
