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
capture mata: mata drop rdistnn_type7()
capture mata: mata drop rdistnn_sample_sd()
capture mata: mata drop rdistnn_scale_sum()
capture mata: mata drop rdistnn_binary_scale()
capture mata: mata drop rdistnn_fit_geometry()
capture mata: mata drop rdistnn_mahal_cov_inv()
capture mata: mata drop rdistnn_mahal_dist()
capture mata: mata drop rdistnn_order()
capture mata: mata drop rdistnn_weights_exact()
capture mata: mata drop rdistnn_weights_asymp()
capture mata: mata drop rdistnn_beta_interval()
capture mata: mata drop rdistnn_ranks()
capture mata: mata drop rdistnn_l2_range_row()
capture mata: mata drop rdistnn_scaled_difference()
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
// L1/L2 distance variants and frozen-rank deletions have independent engine
// selectors; Mahalanobis geometry remains in Mata.
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

// Stable Beta interval probability; use the small upper tail directly above
// the shared lower-endpoint crossover. Preserve every representable mass.
real scalar rdistnn_beta_interval(real scalar a, real scalar b,
    real scalar lower, real scalar upper)
{
    if(lower<=a/(a+b)) return(ibeta(a,b,upper)-ibeta(a,b,lower))
    return(ibetatail(a,b,lower)-ibetatail(a,b,upper))
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
            w_rank[i] = w_rank[i] + rdistnn_beta_interval(l, s-l+1, (i-1)/n, i/n)
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
// R type-7 sample quantile, including even-sample median interpolation.
real scalar rdistnn_type7(real colvector x, real scalar prob)
{
    real scalar h,lo,fraction,left,right
    x=sort(x,1); h=1+(rows(x)-1)*prob; lo=floor(h)
    left=x[lo]; fraction=h-lo
    if(fraction==0 | lo==rows(x)) return(left)
    right=x[lo+1]
    if(left==right) return(left)
    // Match R type 7, including the rounded endpoint products. Subtracting
    // endpoints first can overflow or move a rank-sensitive IQR by one ULP.
    left=(1-fraction)*left
    right=fraction*right
    return(left+right)
}

// Compensated sums are limited to bounded exceptional-scale calculations.
// Products are rounded before entering this shared accumulation convention.
real scalar rdistnn_scale_sum(real colvector x)
{
    real scalar i,total,correction,next_total,value
    total=0; correction=0
    for(i=1;i<=rows(x);i++) {
        value=x[i]; next_total=total+value
        if(abs(total)>=abs(value)) correction=correction+((total-next_total)+value)
        else correction=correction+((value-next_total)+total)
        total=next_total
    }
    return(total+correction)
}

// Mata's ^ uses a transcendental path; construct the binary scale with exact
// multiplications instead. At most eleven iterations cover all double inputs.
real scalar rdistnn_binary_scale(real scalar amplitude)
{
    real scalar exponent,power,factor
    exponent=max((-1074,min((1022,floor(ln(amplitude)/ln(2))))))
    factor=exponent<0 ? .5 : 2
    exponent=abs(exponent); power=1
    while(exponent>0) {
        if(mod(exponent,2)) power=power*factor
        exponent=floor(exponent/2)
        if(exponent>0) factor=factor*factor
    }
    return(power)
}

// Common conservative range gates keep R and Mata on the same arithmetic.
// Binary normalization preserves local differences around large offsets;
// shifting by the first observation keeps the normalized mean well scaled.
real scalar rdistnn_sample_sd(real colvector x)
{
    real scalar amplitude,spread,scale,n
    real colvector z
    n=rows(x); amplitude=max(abs(x))
    if(amplitude<sqrt(4.4942328371557898e307/n)) {
        z=x:-mean(x); spread=max(abs(z))
        if(spread==0) return(0)
        if(spread>=sqrt(2.2250738585072014e-308)*sqrt(n)) {
            return(sqrt(quadcross(z,z)/(n-1)))
        }
    }
    scale=rdistnn_binary_scale(amplitude)
    z=x:/scale; z=z:-z[1]
    z=z:-rdistnn_scale_sum(z)/n
    return(sqrt(rdistnn_scale_sum(z:*z)/(n-1))*scale)
}

// Fit once on the complete training sample. Scaled distances retain scales
// (not inverse scales); Mahalanobis retains its fitted inverse covariance.
real matrix rdistnn_fit_geometry(real matrix X, string scalar distance,
    string scalar rule, real rowvector scalevec, real scalar eps)
{
    real scalar j,n,p,med
    real rowvector scales
    real colvector x
    n=rows(X); p=cols(X)
    if(n<1 | p<1 | hasmissing(X) | missing(eps) | eps<=0) {
        errprintf("Geometry requires finite training data and positive finite scaleeps().\n")
        exit(error(198))
    }
    if(rule=="") rule="sd"
    if(!any(("sd","mad","iqr","range","custom"):==rule)) {
        errprintf("scalerule() must be sd, mad, iqr, range, or custom.\n")
        exit(error(198))
    }
    if((rule=="custom" & (cols(scalevec)!=p | hasmissing(scalevec))) |
        (rule!="custom" & cols(scalevec)>0)) {
        errprintf("scalevec() requires scalerule(custom) and one finite scale per covariate.\n")
        exit(error(198))
    }
    if(distance=="mahalanobis") return(rdistnn_mahal_cov_inv(X))
    if(distance=="l1" | distance=="l2" | distance=="") return(J(0,0,.))
    if(distance!="l1_scaled" & distance!="l2_scaled") {
        errprintf("distance() must be l1, l2, l1_scaled, l2_scaled, or mahalanobis.\n")
        exit(error(198))
    }
    scales=J(1,p,.)
    for(j=1;j<=p;j++) {
        x=X[,j]
        if(rule=="custom") scales[j]=abs(scalevec[j])
        else if(rule=="sd") {
            if(n<2) {
                errprintf("scalerule(sd) requires at least two training observations.\n")
                exit(error(2001))
            }
            scales[j]=rdistnn_sample_sd(x)
        }
        else if(rule=="mad") {
            med=rdistnn_type7(x,.5)
            scales[j]=1.4826*rdistnn_type7(abs(x:-med),.5)
        }
        else if(rule=="iqr") scales[j]=rdistnn_type7(x,.75)-rdistnn_type7(x,.25)
        else scales[j]=max(x)-min(x)
        // Mata's max() skips missing values, so validate before clamping.
        if(missing(scales[j])) {
            errprintf("Fitted scales are nonfinite; rescale covariates.\n")
            exit(error(498))
        }
        scales[j]=max((scales[j],eps))
    }
    if(hasmissing(scales)) {
        errprintf("Fitted scales are nonfinite; rescale covariates.\n")
        exit(error(498))
    }
    return(scales)
}

real colvector rdistnn_ranks(real matrix X, real rowvector x0,
    | string scalar distance, real matrix geometry)
{
    if(args()<3) distance="l2"
    if(args()<4) geometry=rdistnn_fit_geometry(X,distance,"sd",J(1,0,.),1e-12)
    return(rdistnn_order(rdistnn_distance_matrix(X,x0,distance,geometry)))
}

// Shared scalar recovery for scaled L1/L2; valid final differences may fit
// even when the raw subtraction is a Mata missing value.
real scalar rdistnn_scaled_difference(real scalar x,real scalar target,
    real scalar inverse)
{
    real scalar d,half_x,half_t
    d=x-target
    if(missing(d) & abs(inverse)<1) {
        half_x=x*0.5; half_t=target*0.5
        return((half_x-half_t)*(inverse+inverse))
    }
    return(d*inverse)
}

// Recompute an exceptional row with constant accumulator workspace. Ordinary
// rows keep subtract-then-scale, rounded squares and sequential addition.
// R/native use the same 2^1023 intermediate boundary as Mata missing values.
real scalar rdistnn_l2_range_row(real matrix X, real scalar i,
    real rowvector target, real matrix geometry, real scalar scaled)
{
    real scalar c,d,a,scale,ssq,ratio,term
    scale=0; ssq=1
    for(c=1;c<=cols(X);c++) {
        d=X[i,c]-target[c]
        if(scaled) d=rdistnn_scaled_difference(X[i,c],target[c],1/geometry[c])
        a=abs(d)
        if(missing(a)) return(.)
        if(a!=0) {
            if(scale<a) {
                ratio=scale/a
                term=ssq*(ratio*ratio)
                ssq=1+term
                scale=a
            }
            else {
                ratio=a/scale
                term=ratio*ratio
                ssq=ssq+term
            }
        }
    }
    return(scale*sqrt(ssq))
}

// Ordinary distances; frozen geometry can be reused for subsets/new targets.
real matrix rdistnn_distance_matrix(real matrix X, real matrix X0,
    | string scalar distance, real matrix geometry)
{
    real scalar j,m,c,i
    real matrix dmat,A
    real colvector total,exceptional,square,redo
    string scalar engine,plugin_mode
    real scalar native_batch
    st_global("rdistnn_last_distance_engine","")
    engine=st_global("rdistnn_distance_engine")
    if(engine=="") engine="auto"
    if(!any(engine:==("auto","plugin","mata"))) {
        errprintf("rdistnn_distance_engine must be auto, plugin or mata.\n")
        exit(error(198))
    }
    if(rows(X)<1 | cols(X)<1 | cols(X0)!=cols(X) | rows(X0)<1 | hasmissing(X) | hasmissing(X0)) {
        errprintf("Distances require finite aligned covariates and target rows.\n")
        exit(error(198))
    }
    if(args()<3 | distance=="") distance="l2"
    if(args()<4) geometry=rdistnn_fit_geometry(X,distance,"sd",J(1,0,.),1e-12)
    if(distance=="l1_scaled" | distance=="l2_scaled") {
        if(rows(geometry)!=1 | cols(geometry)!=cols(X) | hasmissing(geometry) | any(geometry:<=0)) {
            errprintf("Frozen scale geometry must contain one positive scale per covariate.\n")
            exit(error(498))
        }
    }
    if(distance=="mahalanobis") {
        if(rows(geometry)!=cols(X) | cols(geometry)!=cols(X) | hasmissing(geometry)) {
            errprintf("Frozen Mahalanobis geometry has incompatible dimensions.\n")
            exit(error(498))
        }
    }
    if(any(distance:==("l1","l2","l1_scaled","l2_scaled"))) {
        plugin_mode=st_global("rdistnn_use_plugin")
        // Host matrix transfers dominate small calls. This conservative gate
        // keeps ordinary scalar/small jobs on Mata; explicit requirements
        // always bypass it. The batch kernel itself has no such size limit.
        native_batch=(cols(X)>=16 & rows(X0)>=16 &
            rows(X)*cols(X)*rows(X0)>=2e7)
        if(engine=="plugin" | (engine=="auto" & plugin_mode!="0" &
            (plugin_mode=="1" | native_batch))) {
            dmat=rdistnn_distance_plugin(X,X0,distance,geometry)
            if(rows(dmat)==rows(X) & cols(dmat)==rows(X0)) {
                st_global("rdistnn_last_distance_engine","plugin")
                return(dmat)
            }
            if(engine=="plugin" | plugin_mode=="1") {
                errprintf("Compiled distance API unavailable or failed; install an updated plugin or use automatic/Mata dispatch.\n")
                exit(error(498))
            }
        }
    }
    st_global("rdistnn_last_distance_engine","mata")
    m=rows(X0); dmat=J(rows(X),m,.)
    for(j=1;j<=m;j++) {
        A=X:-X0[j,]
        if(distance=="l1_scaled" | distance=="l2_scaled") A=A:*(1:/geometry)
        // Match the R/native coordinate order, including rounding at near ties.
        // rowsum() may use wider accumulation and therefore change ranks.
        if(distance=="l1" | distance=="l1_scaled") {
            if(distance=="l1_scaled" & hasmissing(A)) {
                // Revisit only failed coordinates; finite ordinary products
                // and their accumulation order are left intact.
                for(c=1;c<=cols(A);c++) {
                    redo=selectindex(A[,c]:>=.)
                    for(i=1;i<=length(redo);i++) {
                        A[redo[i],c]=rdistnn_scaled_difference(
                            X[redo[i],c],X0[j,c],1/geometry[c])
                    }
                }
            }
            total=J(rows(X),1,0)
            for(c=1;c<=cols(A);c++) total=total:+abs(A[,c])
            dmat[,j]=total
        }
        else if(distance=="l2" | distance=="l2_scaled") {
            total=J(rows(X),1,0)
            exceptional=J(rows(X),1,0)
            for(c=1;c<=cols(A);c++) {
                square=A[,c]:*A[,c]
                // Missing squares propagate to total; form row masks only
                // when a square might have underflowed, matching the R gate.
                if(any(square:<2.2250738585072014e-308)) {
                    exceptional=exceptional :|
                        ((A[,c]:!=0) :& (square:<2.2250738585072014e-308))
                }
                total=total:+square
            }
            dmat[,j]=sqrt(total)
            redo=selectindex(exceptional :| (total:>=.))
            // selectindex() returns a 1 x 0 empty vector for scalar false.
            for(i=1;i<=length(redo);i++) {
                dmat[redo[i],j]=rdistnn_l2_range_row(X,redo[i],X0[j,],
                    geometry,distance=="l2_scaled")
            }
        }
        else if(distance=="mahalanobis") dmat[,j]=rdistnn_mahal_dist(X,X0[j,],geometry)
        else {
            errprintf("Unknown distance '%s'.\n",distance)
            exit(error(198))
        }
    }
    if(hasmissing(dmat)) {
        errprintf("Distances overflowed; rescale covariates or increase scaleeps().\n")
        exit(error(498))
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
                                  | string scalar distance, real matrix geometry)
{
    if (args() < 3) distance = "l2"
    if (args() < 4) geometry=rdistnn_fit_geometry(X,distance,"sd",J(1,0,.),1e-12)
    return(rdistnn_ranks_from_dmat(
        rdistnn_distance_matrix(X, X0, distance, geometry)))
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

    if (args() < 7 | rows(Sinv) == 0) Sinv=rdistnn_fit_geometry(X,distance,"sd",J(1,0,.),1e-12)
    dist=rdistnn_distance_matrix(X,x0,distance,Sinv)
    ranks = rdistnn_order(dist)
    return(rdistnn_weights_from_ranks(ranks, k, s, method))
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

    loo = rdistnn_loo_from_ranks(ranks, signal, k, s, method, 1)
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

    loo = rdistnn_loo_from_ranks(ranks, signal, k, s, method, 1)
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

// Shared optional-native distance bridge, separate from geometry fitting.
findfile "_rdistnn_distance_plugin.ado"
quietly run "`r(fn)'"

// Keep the shared arbitrary-table engine separate from geometry and weights.
findfile "_rdistnn_rank_loo.ado"
quietly run "`r(fn)'"
