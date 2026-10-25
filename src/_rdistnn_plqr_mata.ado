*! version 0.1.5  2026-10-25
version 17
// Adaptive PLQR predicate: retain certified ordinary comparisons, otherwise
// sum the binary64 products exactly with separate base-2^26 magnitudes.
foreach fn in indicator exact_indicator exact_product exact_piece {
    capture mata: mata drop rdistnn_plqr_`fn'()
}
mata:
void rdistnn_plqr_exact_piece(real colvector limbs, real scalar piece,
    real scalar bit)
{
    real scalar index,value,next_value,carry
    if(piece==0) return
    if(bit<0) rdistnn_fail(498,"Internal PLQR comparison exponent is invalid; report a reproducible example.")
    index=floor(bit/26)+1
    value=piece*rdistnn_power2(mod(bit,26))
    while(value!=0) {
        if(index>rows(limbs)) {
            rdistnn_fail(498,"Internal PLQR exact comparison capacity exceeded.")
        }
        next_value=limbs[index]+value
        carry=floor(next_value/67108864)
        limbs[index]=next_value-carry*67108864
        value=carry; index=index+1
    }
}

void rdistnn_plqr_exact_product(real colvector positive,
    real colvector negative, real scalar a, real scalar b)
{
    real rowvector left,right,xs,ys
    real scalar x,y,xnext,ynext,i,j,product,high,low,start,bit
    if(a==0 | b==0) return
    left=rdistnn_binary_parts(a); right=rdistnn_binary_parts(b)
    x=abs(left[1])*4503599627370496
    y=abs(right[1])*4503599627370496
    xs=ys=J(1,3,0)
    for(i=1;i<=3;i++) {
        xnext=floor(x/67108864); ynext=floor(y/67108864)
        xs[i]=x-xnext*67108864; ys[i]=y-ynext*67108864
        x=xnext; y=ynext
    }
    // Origin -2252 covers two minimum operand powers, with the two
    // 52-bit significand shifts. Split products before the bit shift so
    // every multiplication/addition is an exact integer below 2^52.
    start=left[2]+right[2]-104-(-2252)
    for(i=1;i<=3;i++) for(j=1;j<=3;j++) {
        product=xs[i]*ys[j]
        high=floor(product/67108864); low=product-high*67108864
        bit=start+26*(i+j-2)
        if((a<0)==(b<0)) {
            rdistnn_plqr_exact_piece(positive,low,bit)
            rdistnn_plqr_exact_piece(positive,high,bit+26)
        }
        else {
            rdistnn_plqr_exact_piece(negative,low,bit)
            rdistnn_plqr_exact_piece(negative,high,bit+26)
        }
    }
}

real scalar rdistnn_plqr_exact_indicator(real rowvector D,
    real colvector theta, real scalar g, real scalar y)
{
    real colvector positive,negative
    real scalar j
    // 170 limbs cover [-2252,2168). At most INT_MAX+2 products, each
    // smaller than 2^2048, have magnitude below 2^2080 in either sign.
    positive=negative=J(170,1,0)
    for(j=1;j<=cols(D);j++) {
        rdistnn_plqr_exact_product(positive,negative,D[j],theta[j])
    }
    rdistnn_plqr_exact_product(positive,negative,g,1)
    // Keep y and -1 separate: negating a finite negative extreme can
    // exceed Stata's positive ceiling before the exact product is formed.
    rdistnn_plqr_exact_product(positive,negative,y,-1)
    for(j=170;j>=1;j--) {
        if(positive[j]>negative[j]) return(1)
        if(positive[j]<negative[j]) return(0)
    }
    return(1)
}

real colvector rdistnn_plqr_indicator(real matrix D, real colvector theta,
    real colvector g, real colvector y)
{
    real scalar n,q,i,j,risk,product,margin,bound,gate,tiny
    real colvector threshold,below
    n=rows(D); q=cols(D)
    if(n<1 | q<1 | q>2147483647 | rows(theta)!=q | rows(g)!=n |
        rows(y)!=n | hasmissing(D) | hasmissing(theta) |
        hasmissing(g) | hasmissing(y)) {
        rdistnn_fail(198,"PLQR comparisons require finite aligned regressors, coefficients, nuisance values and outcomes.")
    }
    gate=rdistnn_power2(1022); tiny=rdistnn_power2(-1074)
    threshold=D*theta:+g
    below=(y:<=threshold)
    for(i=1;i<=n;i++) {
        // Match the R/native risk bound and materialized product order.
        risk=abs(g[i])+abs(y[i])
        for(j=1;j<=q;j++) {
            product=D[i,j]*theta[j]
            risk=risk+abs(product)
        }
        if(risk<gate & threshold[i]<.) {
            margin=threshold[i]-y[i]
            bound=8*(q+3)*epsilon(1)*risk+32*(q+2)*tiny
            if(margin<. & abs(margin)>bound) continue
        }
        below[i]=rdistnn_plqr_exact_indicator(D[i,],theta,g[i],y[i])
    }
    return(below)
}
end

// Fresh-command IO with function-local scratch.
capture mata: mata drop rdistnn_plqr_fit_ado()
mata:
void rdistnn_plqr_fit_ado()
{
    real scalar ado_k,ado_s,ado_scaleeps,ado_tau,_pl_neff
    string scalar ado_V,ado_b,ado_diagnostics,ado_distance,ado_dvars,ado_gvar,ado_mdvars,ado_method
    string scalar ado_model_distances,ado_model_geometry,ado_model_groups,ado_model_ranks
    string scalar ado_model_signals,ado_outcome,ado_scalerule,ado_scalevec,ado_targets,ado_theta_mat
    string scalar ado_thetaeval,ado_touse,ado_x0,ado_xvars
    real matrix _pl_D,_pl_DX,_pl_V,_pl_X,_pl_below,_pl_g0,_pl_geometry,_pl_groups,_pl_md,_pl_moment,_pl_ranks,_pl_resid,_pl_score,_pl_theta,_pl_w,_pl_x0,_pl_y
    ado_k=strtoreal(st_local("k"))
    ado_s=strtoreal(st_local("s"))
    ado_scaleeps=strtoreal(st_local("scaleeps"))
    ado_tau=strtoreal(st_local("tau"))
    ado_V=st_local("V")
    ado_b=st_local("b")
    ado_diagnostics=st_local("diagnostics")
    ado_distance=st_local("distance")
    ado_dvars=st_local("dvars")
    ado_gvar=st_local("gvar")
    ado_mdvars=st_local("mdvars")
    ado_method=st_local("method")
    ado_model_distances=st_local("model_distances")
    ado_model_geometry=st_local("model_geometry")
    ado_model_groups=st_local("model_groups")
    ado_model_ranks=st_local("model_ranks")
    ado_model_signals=st_local("model_signals")
    ado_outcome=st_local("outcome")
    ado_scalerule=st_local("scalerule")
    ado_scalevec=st_local("scalevec")
    ado_targets=st_local("targets")
    ado_theta_mat=st_local("theta_mat")
    ado_thetaeval=st_local("thetaeval")
    ado_touse=st_local("touse")
    ado_x0=st_local("x0")
    ado_xvars=st_local("xvars")

    _pl_X      = st_data(., ado_xvars, ado_touse)
    _pl_y      = st_data(., ado_outcome, ado_touse)
    _pl_D      = st_data(., ado_dvars, ado_touse)
    _pl_g0     = st_data(., ado_gvar, ado_touse)
    _pl_md     = st_data(., ado_mdvars, ado_touse)
    _pl_x0     = strtoreal(tokens(ado_x0))
    _pl_theta  = strtoreal(tokens(ado_thetaeval))'
    // Match the R reference and the other estimators: rank the reported
    // distance matrix, with original observation IDs ordering exact ties.
    _pl_geometry = rdistnn_fit_geometry(_pl_X,ado_distance,ado_scalerule,strtoreal(tokens(ado_scalevec)),ado_scaleeps)
    _pl_DX     = rdistnn_distance_matrix(_pl_X, _pl_x0, ado_distance, _pl_geometry)
    _pl_ranks  = rdistnn_ranks_from_dmat(_pl_DX)
    _pl_groups=rdistnn_distance_groups(_pl_ranks,_pl_DX,1)
    _pl_w      = rdistnn_weights_from_rmat(_pl_ranks, ado_k, ado_s, ado_method,_pl_DX,_pl_groups)
    _pl_below = rdistnn_plqr_indicator(_pl_D,_pl_theta,_pl_g0,_pl_y)
    _pl_resid  = ado_tau :- _pl_below
    _pl_score  = (_pl_D :- _pl_md) :* (_pl_resid * J(1, cols(_pl_D), 1))
    _pl_moment = rdistnn_weighted_signals(_pl_w, _pl_score)
    _pl_V      = rdistnn_jk_cov_from_ranks(_pl_ranks, _pl_score, ado_k, ado_s, ado_method,_pl_DX,_pl_groups)
    _pl_neff   = 1 / sum(_pl_w :^2)
    st_matrix(ado_b, _pl_moment)
    st_matrix(ado_V, _pl_V)
    st_matrix(ado_theta_mat, _pl_theta')
    st_numscalar("r(neff)", _pl_neff)
    st_matrix(ado_targets,_pl_x0)
    st_matrix(ado_model_groups,_pl_groups); st_matrix(ado_model_signals,_pl_score); st_matrix(ado_model_ranks,_pl_ranks); st_matrix(ado_model_distances,_pl_DX)
    if(any(ado_distance:==("mahalanobis","l1_scaled","l2_scaled"))) st_matrix(ado_model_geometry,_pl_geometry)
    st_matrix(ado_diagnostics,rdistnn_diagnostics("plqr_dnn",_pl_X,_pl_x0,_pl_score,_pl_w,_pl_DX,ado_k,ado_s,_pl_moment,J(0,1,.),J(0,1,.)))

}
end
