*! version 0.1.5  2026-10-25
version 17
// Signed two-scale rank tables and generic, frozen-rank deletion inference.
foreach fn in coefficients table weights cov diagnostics {
    capture mata: mata drop rdistnn_tdnn_`fn'()
}
mata:
// Shared coefficient arithmetic avoids cancellation at nearby scales or
// large bias dimension. Validation remains with the table's public contract.
real rowvector rdistnn_tdnn_coefficients(real scalar s1, real scalar s2,
    real scalar dim)
{
    real scalar w1
    w1=-1/expm1((2/dim)*ln1p((s2-s1)/s1))
    return((w1,1-w1))
}

real colvector rdistnn_tdnn_table(real scalar n, real scalar s1,
    real scalar s2, real scalar dim, string scalar method, real scalar normalize)
{
    real scalar r,total
    real rowvector coefficients
    real colvector a,b,w
    if(hasmissing((n,s1,s2,dim,normalize)) | n<1 | n!=floor(n) |
        s1<1 | s1!=floor(s1) | s2<=s1 | s2!=floor(s2) | s2>n |
        dim<1 | dim!=floor(dim) | (normalize!=0 & normalize!=1) |
        (method!="exact" & method!="asymp")) {
        rdistnn_fail(198,"TDNN requires integer 1 <= s1 < s2 <= n, positive dim, exact/asymp and logical normalization.")
    }
    a=b=J(n,1,0)
    if(method=="exact") {
        // The same portable exact components used by DNN and R TDNN retain
        // signed cancellation consistently before the two-scale combination.
        a=rdistnn_exact_table(n,1,s1)
        b=rdistnn_exact_table(n,1,s2)
    }
    else {
        if(s1==1) a=J(n,1,1/n)
        for(r=1;r<=n;r++) {
            if(s1>1) a[r]=rdistnn_beta_interval(1,s1,(r-1)/n,r/n)
            b[r]=rdistnn_beta_interval(1,s2,(r-1)/n,r/n)
        }
    }
    coefficients=rdistnn_tdnn_coefficients(s1,s2,dim)
    w=coefficients[1]*a+coefficients[2]*b
    if(hasmissing(w)) {
        rdistnn_fail(498,"TDNN scale combination produced nonfinite weights; choose less extreme localization orders.")
    }
    if(method=="exact" | normalize) {
        total=rdistnn_sum(w)
        if(missing(total) | abs(total)<100*epsilon(1)) rdistnn_fail(498,"TDNN normalization is unavailable; choose less extreme localization orders.")
        if(abs(total-1)>1e-10) w=w/total
    }
    return(w)
}

real matrix rdistnn_tdnn_weights(real matrix ranks, real scalar s1,
    real scalar s2, real scalar dim, string scalar method, real scalar normalize,
    | real matrix DX, real matrix prepared_groups)
{
    real colvector table
    real matrix W
    real scalar j,n
    if(args()<8) prepared_groups=J(0,0,.)
    n=rows(ranks)
    table=rdistnn_tdnn_table(n,s1,s2,dim,method,normalize)
    if(cols(ranks)<1 | hasmissing(ranks)) rdistnn_fail(198,"TDNN requires at least one target column with finite ranks.")
    W=J(n,cols(ranks),.)
    for(j=1;j<=cols(ranks);j++) {
        if(any(sort(ranks[,j],1):!=(1::n))) {
            rdistnn_fail(198,"TDNN ranks must be permutations of 1..n.")
        }
        W[,j]=table[ranks[,j]]
    }
    if(args()>=7 & rows(DX)>0) W=rdistnn_distance_weights(ranks,W,DX,prepared_groups)
    return(W)
}

real matrix rdistnn_tdnn_cov(real matrix ranks, real matrix signal,
    real scalar s1, real scalar s2, real scalar dim, string scalar method,
    real scalar normalize, | real matrix DX, real matrix prepared_groups)
{
    real scalar n
    real matrix loo
    if(args()<9) prepared_groups=J(0,0,.)
    n=rows(ranks)
    if(args()>=8 & rows(DX)>0) {
        loo=rdistnn_loo_from_table(ranks,signal,
            rdistnn_tdnn_table(n-1,s1,s2,dim,method,normalize),1,"",DX,prepared_groups)
    }
    else loo=rdistnn_loo_from_table(ranks,signal,
        rdistnn_tdnn_table(n-1,s1,s2,dim,method,normalize),1)
    return(rdistnn_covariance(loo,(n-1)/n,1))
}

// No signed-probability quantiles, variances or effective sample sizes.
real matrix rdistnn_tdnn_diagnostics(real matrix W, real matrix DX)
{
    real scalar j,total
    real colvector w
    real matrix out
    if(rows(W)!=rows(DX) | cols(W)!=cols(DX)) rdistnn_fail(503,"TDNN weights and distances must have the same dimensions.")
    out=J(cols(W),8,.)
    for(j=1;j<=cols(W);j++) {
        w=W[,j]; total=quadsum(abs(w))
        out[j,]=(quadsum(w),quadsum(select(w,w:>0)),
            -quadsum(select(w,w:<0)),total,1/quadsum(w:^2),max(abs(w)),
            sum(w:!=0),cross(abs(w),DX[,j])/total)
    }
    return(out)
}
end

// Fresh-command IO with function-local scratch.
capture mata: mata drop rdistnn_tdnn_weights_ado()
capture mata: mata drop rdistnn_tdnn_signal_ado()
mata:
void rdistnn_tdnn_weights_ado()
{
    real scalar ado_dim,ado_n,ado_normalize,ado_s1,ado_s2,ado_scaleeps
    string scalar ado_DX,ado_R,ado_W,ado_diag,ado_distance,ado_geom,ado_groups,ado_method,ado_row_ids
    string scalar ado_scalerule,ado_scalevec,ado_table,ado_targets,ado_touse,ado_varlist
    real matrix _td_DX,_td_R,_td_W,_td_X,_td_X0,_td_geom,_td_groups,_td_table
    ado_dim=strtoreal(st_local("dim"))
    ado_n=strtoreal(st_local("n"))
    ado_normalize=strtoreal(st_local("normalize"))
    ado_s1=strtoreal(st_local("s1"))
    ado_s2=strtoreal(st_local("s2"))
    ado_scaleeps=strtoreal(st_local("scaleeps"))
    ado_DX=st_local("DX")
    ado_R=st_local("R")
    ado_W=st_local("W")
    ado_diag=st_local("diag")
    ado_distance=st_local("distance")
    ado_geom=st_local("geom")
    ado_groups=st_local("groups")
    ado_method=st_local("method")
    ado_row_ids=st_local("row_ids")
    ado_scalerule=st_local("scalerule")
    ado_scalevec=st_local("scalevec")
    ado_table=st_local("table")
    ado_targets=st_local("targets")
    ado_touse=st_local("touse")
    ado_varlist=st_local("varlist")

    _td_X=st_data(.,tokens(ado_varlist),ado_touse)
    _td_X0=st_matrix(ado_targets)
    if(hasmissing(_td_X0)) rdistnn_fail(198,"x0matrix() must contain only finite target values; supply one target per row.")
    _td_geom=rdistnn_fit_geometry(_td_X,ado_distance,ado_scalerule,strtoreal(tokens(ado_scalevec)),ado_scaleeps)
    _td_DX=rdistnn_distance_matrix(_td_X,_td_X0,ado_distance,_td_geom)
    _td_R=rdistnn_ranks_from_dmat(_td_DX)
    _td_groups=rdistnn_distance_groups(_td_R,_td_DX,1)
    rdistnn_distance_warning(_td_groups)
    _td_table=rdistnn_tdnn_table(ado_n,ado_s1,ado_s2,ado_dim,ado_method,ado_normalize)
    _td_W=rdistnn_tdnn_weights(_td_R,ado_s1,ado_s2,ado_dim,ado_method,ado_normalize,_td_DX,_td_groups)
    st_matrix(ado_geom,_td_geom)
    st_matrix(ado_row_ids,select((1::st_nobs()),st_data(.,ado_touse)))
    st_matrix(ado_W,_td_W); st_matrix(ado_R,_td_R); st_matrix(ado_DX,_td_DX)
    st_matrix(ado_groups,_td_groups)
    st_matrix(ado_table,_td_table)
    st_matrix(ado_diag,rdistnn_tdnn_diagnostics(_td_W,_td_DX))
}

void rdistnn_tdnn_signal_ado()
{
    string scalar ado_W,ado_b,ado_outcome,ado_signals,ado_touse
    real matrix _td_y
    ado_W=st_local("W")
    ado_b=st_local("b")
    ado_outcome=st_local("outcome")
    ado_signals=st_local("signals")
    ado_touse=st_local("touse")

    _td_y=st_data(.,ado_outcome,ado_touse)
    st_matrix(ado_signals,_td_y)
    st_matrix(ado_b,rdistnn_weighted_signals(_td_y,st_matrix(ado_W)))
}
end
