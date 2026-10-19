version 17
// Fitted geometry, distance evaluation and stable rank ordering.
capture mata: mata drop rdistnn_l2dist()
capture mata: mata drop rdistnn_mahal_cov_inv()
capture mata: mata drop rdistnn_mahal_dist()
capture mata: mata drop rdistnn_order()
capture mata: mata drop rdistnn_type7()
capture mata: mata drop rdistnn_sample_sd()
capture mata: mata drop rdistnn_fit_geometry()
capture mata: mata drop rdistnn_ranks()
capture mata: mata drop rdistnn_scaled_difference()
capture mata: mata drop rdistnn_l2_range_row()
capture mata: mata drop rdistnn_distance_matrix()
capture mata: mata drop rdistnn_ranks_from_dmat()
capture mata: mata drop rdistnn_ranks_matrix()
capture mata: mata drop rdistnn_physical_distances()
mata:
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

// Common conservative range gates keep R and Mata on the same arithmetic.
// Binary normalization preserves local differences around large offsets;
// shifting by the first observation keeps the normalized mean well scaled.
real scalar rdistnn_sample_sd(real colvector x)
{
    real scalar amplitude,spread,scale,n
    real colvector z
    n=rows(x); amplitude=max(abs(x))
    if(amplitude<sqrt(4.4942328371557898e307/n)) {
        z=x:-rdistnn_signal_level(x); z=z:-mean(z); spread=max(abs(z))
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
// (not inverse scales); Mahalanobis retains a centered-design whitening factor.
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
    if(distance=="mahalanobis") return(rdistnn_mahal_geometry(X))
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
        if(rows(geometry)==cols(X)+2 & cols(geometry)==cols(X) & !hasmissing(geometry)) {
            if(geometry[1,1]<=0) {
                errprintf("Frozen Mahalanobis geometry requires a positive scale.\n")
                exit(error(498))
            }
            return(rdistnn_distance_matrix(rdistnn_mahal_project(X,geometry),
                rdistnn_mahal_project(X0,geometry),"l2",J(0,0,.)))
        }
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

// Physical radius always uses ordinary original-coordinate Euclidean distance.
// Its internal dispatch must not overwrite the fitted metric's engine marker.
real colvector rdistnn_physical_distances(real matrix X, real rowvector x0)
{
    string scalar engine
    real colvector distances
    engine=st_global("rdistnn_last_distance_engine")
    distances=rdistnn_distance_matrix(X,x0,"l2",J(0,0,.))
    st_global("rdistnn_last_distance_engine",engine)
    return(distances)
}
end
