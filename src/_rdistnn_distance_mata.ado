*! version 0.1.5  2026-10-25
version 17
// Fitted geometry, distance evaluation and stable rank ordering.
capture mata: mata drop rdistnn_order()
capture mata: mata drop rdistnn_type7()
capture mata: mata drop rdistnn_sample_sd()
capture mata: mata drop rdistnn_fit_geometry()
capture mata: mata drop rdistnn_scaled_difference()
capture mata: mata drop rdistnn_distance_power2()
capture mata: mata drop rdistnn_scaled_norm_row()
capture mata: mata drop rdistnn_l2_range_row()
capture mata: mata drop rdistnn_distance_matrix()
capture mata: mata drop rdistnn_ranks_from_dmat()
capture mata: mata drop rdistnn_physical_distances()
mata:


// ── rdistnn_order ─────────────────────────────────────────────────────
// Return the integer rank order (1 = nearest) of each observation
// by distance to x0, breaking exact represented ties by original row order.
//
// Returns: n-vector of integer ranks.
//
real colvector rdistnn_order(real colvector dist2)
{
    real colvector  idx, rnk
    real scalar     n, i

    n   = rows(dist2)
    idx = order((dist2, (1::n)), (1, 2))
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
    if(!hasmissing(abs(x)) & amplitude<sqrt(4.4942328371557898e307/n)) {
        z=x:-rdistnn_signal_level(x); z=z:-mean(z); spread=max(abs(z))
        if(spread==0) return(0)
        if(spread>=sqrt(2.2250738585072014e-308)*sqrt(n)) {
            return(sqrt(quadcross(z,z)/(n-1)))
        }
    }
    // Computed negative inputs can exceed Stata's positive finite ceiling;
    // abs() then looks missing, so use a bounded exact binary scale directly.
    if(hasmissing(abs(x))) scale=rdistnn_power2(1022)
    else scale=rdistnn_binary_scale(amplitude)
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
        rdistnn_fail(198,"Geometry requires finite training data and positive finite scaleeps().")
    }
    if(rule=="") rule="sd"
    if(!any(("sd","mad","iqr","range","custom"):==rule)) {
        rdistnn_fail(198,"scalerule() must be sd, mad, iqr, range, or custom.")
    }
    if((rule=="custom" & (cols(scalevec)!=p | hasmissing(scalevec))) |
        (rule!="custom" & cols(scalevec)>0)) {
        rdistnn_fail(198,"scalevec() requires scalerule(custom) and one finite scale per covariate.")
    }
    if(distance=="mahalanobis") return(rdistnn_mahal_geometry(X))
    if(distance=="l1" | distance=="l2" | distance=="") return(J(0,0,.))
    if(distance!="l1_scaled" & distance!="l2_scaled") {
        rdistnn_fail(198,"distance() must be l1, l2, l1_scaled, l2_scaled, or mahalanobis.")
    }
    scales=J(1,p,.)
    for(j=1;j<=p;j++) {
        x=X[,j]
        if(rule=="custom") scales[j]=abs(scalevec[j])
        else if(rule=="sd") {
            if(n<2) {
                rdistnn_fail(2001,"scalerule(sd) requires at least two training observations.")
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
            rdistnn_fail(498,"Fitted scales are nonfinite; rescale covariates.")
        }
        scales[j]=max((scales[j],eps))
    }
    if(hasmissing(scales)) {
        rdistnn_fail(498,"Fitted scales are nonfinite; rescale covariates.")
    }
    return(scales)
}


// Shared scalar recovery for scaled L1/L2; valid final differences may fit
// even when the raw subtraction is a Mata missing value.
real scalar rdistnn_scaled_difference(real scalar x,real scalar target,
    real scalar inverse)
{
    real scalar d,half_x,half_t,half
    d=x-target
    if(abs(d)>=. & abs(inverse)<1) {
        half_x=x*0.5; half_t=target*0.5
        half=half_x-half_t
        if(abs(half)>=.) return((x*.25-target*.25)*(4*inverse))
        return(half*(inverse+inverse))
    }
    return(d*inverse)
}

// Split the restoration so a representable final subnormal is rounded once.
// Match R truncation and C++ integer division, including negative powers.
real scalar rdistnn_distance_power2(real scalar value, real scalar power)
{
    real scalar half
    if(value==0) return(value)
    half=trunc(power/2)
    return((value*rdistnn_power2(half))*rdistnn_power2(power-half))
}

// Exceptional scaled products retain binary exponents through the complete
// row norm. Ordinary rows continue using their original rounded products.
real scalar rdistnn_scaled_norm_row(real matrix X, real scalar i,
    real rowvector target, real rowvector inverse, real scalar l2,
    | real scalar squared)
{
    real scalar c,raw,shift,m,e,scale_m,scale_e,total,ratio,square,term
    real scalar active,rounded_max,rounded,norm
    real rowvector left,right
    if(args()<6) squared=0
    scale_m=total=0; scale_e=0; active=0; rounded_max=0
    for(c=1;c<=cols(X);c++) {
        raw=X[i,c]-target[c]; shift=0
        if(abs(raw)>=.) {
            raw=X[i,c]*.25-target[c]*.25
            shift=2
        }
        if(raw==0 | inverse[c]==0) continue
        if(missing(raw) | missing(inverse[c])) return(abs(raw*inverse[c]))
        active=active+1
        rounded=abs(rdistnn_scaled_difference(X[i,c],target[c],inverse[c]))
        // Mata max() ignores missing, whereas a coordinate overflow must
        // remain an explicit distance overflow after norm recovery.
        if(missing(rounded)) return(rounded)
        rounded_max=max((rounded_max,rounded))
        left=rdistnn_binary_parts(abs(raw))
        right=rdistnn_binary_parts(abs(inverse[c]))
        while(left[1]>=1) {
            left[1]=left[1]/2; left[2]=left[2]+1
        }
        while(right[1]>=1) {
            right[1]=right[1]/2; right[2]=right[2]+1
        }
        m=left[1]*right[1]; e=left[2]+right[2]+shift
        if(m<.5) {
            m=m*2; e=e-1
        }
        if(scale_m==0) {
            scale_m=m; scale_e=e; total=1
        }
        else if(e>scale_e | (e==scale_e & m>scale_m)) {
            ratio=rdistnn_distance_power2(scale_m/m,scale_e-e)
            square=ratio*ratio
            term=total*(l2 ? square : ratio)
            total=1+term
            scale_m=m; scale_e=e
        }
        else {
            ratio=rdistnn_distance_power2(m/scale_m,e-scale_e)
            square=ratio*ratio
            total=total+(l2 ? square : ratio)
        }
    }
    if(scale_m==0) return(0)
    // Preserve direct product rounding for a single active coordinate,
    // including rows padded with zeros at subnormal rounding boundaries.
    if(active==1) return(squared ? rounded_max*rounded_max : rounded_max)
    if(squared) return(rdistnn_distance_power2((scale_m*total)*scale_m,2*scale_e))
    norm=rdistnn_distance_power2(scale_m*(l2 ? sqrt(total) : total),scale_e)
    if(missing(norm)) return(norm)
    // Monotonic rounding preserves the exact norm's coordinate lower bound.
    return(max((norm,rounded_max)))
}

// Recompute an exceptional row with constant accumulator workspace. Ordinary
// rows keep subtract-then-scale, rounded squares and sequential addition.
// R/native use the same 2^1023 intermediate boundary as Mata missing values.
real scalar rdistnn_l2_range_row(real matrix X, real scalar i,
    real rowvector target, real matrix geometry, real scalar scaled)
{
    real scalar c,d,a,scale,ssq,ratio,term
    if(scaled) return(rdistnn_scaled_norm_row(X,i,target,1:/geometry,1))
    scale=0; ssq=1
    for(c=1;c<=cols(X);c++) {
        d=X[i,c]-target[c]
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
    real scalar j,m,c,i,inverse
    real matrix dmat,A
    real colvector total,exceptional,square,redo,raw
    string scalar engine,plugin_mode
    real scalar native_batch
    st_global("rdistnn_last_distance_engine","")
    engine=st_global("rdistnn_distance_engine")
    if(engine=="") engine="auto"
    if(!any(engine:==("auto","plugin","mata"))) {
        rdistnn_fail(198,"rdistnn_distance_engine must be auto, plugin or mata.")
    }
    if(rows(X)<1 | cols(X)<1 | cols(X0)!=cols(X) | rows(X0)<1 | hasmissing(X) | hasmissing(X0)) {
        rdistnn_fail(198,"Distances require finite aligned covariates and target rows.")
    }
    if(args()<3 | distance=="") distance="l2"
    if(args()<4) geometry=rdistnn_fit_geometry(X,distance,"sd",J(1,0,.),1e-12)
    if(distance=="l1_scaled" | distance=="l2_scaled") {
        if(rows(geometry)!=1 | cols(geometry)!=cols(X) | hasmissing(geometry) | any(geometry:<=0)) {
            rdistnn_fail(498,"Frozen scale geometry must contain one positive scale per covariate; refit the model.")
        }
    }
    if(distance=="mahalanobis") {
        if(rows(geometry)!=cols(X)+2 | cols(geometry)!=cols(X) |
            hasmissing(geometry) | geometry[1,1]<=0) {
            rdistnn_fail(498,"Mahalanobis distance requires current fitted whitening geometry; refit the model.")
        }
        return(rdistnn_distance_matrix(rdistnn_mahal_project(X,geometry),
            rdistnn_mahal_project(X0,geometry),"l2",J(0,0,.)))
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
                rdistnn_fail(498,"Compiled distance API is unavailable; install the current plugin or set rdistnn_distance_engine mata.")
            }
        }
    }
    st_global("rdistnn_last_distance_engine","mata")
    m=rows(X0); dmat=J(rows(X),m,.)
    for(j=1;j<=m;j++) {
        A=X:-X0[j,]
        exceptional=J(rows(X),1,0)
        if(distance=="l1_scaled" | distance=="l2_scaled") {
            for(c=1;c<=cols(A);c++) {
                raw=A[,c]; inverse=1/geometry[c]
                exceptional=exceptional :| (abs(raw):>=.)
                A[,c]=raw:*inverse
                exceptional=exceptional :|
                    ((raw:!=0) :& (inverse!=0) :& (abs(A[,c]):<2.2250738585072014e-308))
            }
        }
        // Match the R/native coordinate order, including rounding at near ties.
        // rowsum() may use wider accumulation and therefore change ranks.
        if(distance=="l1" | distance=="l1_scaled") {
            total=J(rows(X),1,0)
            for(c=1;c<=cols(A);c++) total=total:+abs(A[,c])
            dmat[,j]=total
            if(distance=="l1_scaled") {
                redo=selectindex(exceptional)
                for(i=1;i<=length(redo);i++) {
                    dmat[redo[i],j]=rdistnn_scaled_norm_row(X,redo[i],X0[j,],1:/geometry,0)
                }
            }
        }
        else if(distance=="l2" | distance=="l2_scaled") {
            total=J(rows(X),1,0)
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
        else {
            rdistnn_fail(198,"Unknown distance(): "+rdistnn_message_value(distance)+". Use l1, l2, l1_scaled, l2_scaled, or mahalanobis.")
        }
    }
    if(hasmissing(dmat)) {
        rdistnn_fail(498,"Distances overflowed; rescale covariates or increase scaleeps().")
    }
    return(dmat)
}

// Convert each distance column to a full rank permutation.
real matrix rdistnn_ranks_from_dmat(real matrix dmat)
{
    real scalar j
    real matrix ranks

    if (rows(dmat) < 1 | cols(dmat) < 1) {
        rdistnn_fail(198,"rdistnn_ranks_from_dmat: dmat must be nonempty.")
    }
    ranks = J(rows(dmat), cols(dmat), .)
    for (j = 1; j <= cols(dmat); j++) {
        ranks[,j] = rdistnn_order(dmat[,j])
    }
    return(ranks)
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
