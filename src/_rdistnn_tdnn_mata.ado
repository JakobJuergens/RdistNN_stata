*! version 0.1.0  2026-10-18
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
        errprintf("TDNN requires integer 1 <= s1 < s2 <= n, positive dim, exact/asymp and logical normalization.\n")
        exit(error(198))
    }
    a=b=J(n,1,0)
    if(method=="exact") {
        // The same portable exact components used by DNN and R TDNN retain
        // signed cancellation consistently before the two-scale combination.
        a=rdistnn_exact_table(n,1,s1)
        b=rdistnn_exact_table(n,1,s2)
    }
    else {
        for(r=1;r<=n;r++) {
            a[r]=rdistnn_beta_interval(1,s1,(r-1)/n,r/n)
            b[r]=rdistnn_beta_interval(1,s2,(r-1)/n,r/n)
        }
    }
    coefficients=rdistnn_tdnn_coefficients(s1,s2,dim)
    w=coefficients[1]*a+coefficients[2]*b
    if(hasmissing(w)) {
        errprintf("TDNN scale combination produced nonfinite weights.\n")
        exit(error(498))
    }
    if(method=="exact" | normalize) {
        total=rdistnn_sum(w)
        if(missing(total) | abs(total)<100*epsilon(1)) exit(error(498))
        if(abs(total-1)>1e-10) w=w/total
    }
    return(w)
}

real matrix rdistnn_tdnn_weights(real matrix ranks, real scalar s1,
    real scalar s2, real scalar dim, string scalar method, real scalar normalize)
{
    real colvector table
    real matrix W
    real scalar j,n
    n=rows(ranks)
    table=rdistnn_tdnn_table(n,s1,s2,dim,method,normalize)
    if(cols(ranks)<1 | hasmissing(ranks)) exit(error(198))
    W=J(n,cols(ranks),.)
    for(j=1;j<=cols(ranks);j++) {
        if(any(sort(ranks[,j],1):!=(1::n))) {
            errprintf("TDNN ranks must be permutations of 1..n.\n")
            exit(error(198))
        }
        W[,j]=table[ranks[,j]]
    }
    return(W)
}

real matrix rdistnn_tdnn_cov(real matrix ranks, real matrix signal,
    real scalar s1, real scalar s2, real scalar dim, string scalar method,
    real scalar normalize)
{
    real scalar n
    real matrix loo
    n=rows(ranks)
    loo=rdistnn_loo_from_table(ranks,signal,
        rdistnn_tdnn_table(n-1,s1,s2,dim,method,normalize),1)
    return(rdistnn_covariance(loo,(n-1)/n,1))
}

// No signed-probability quantiles, variances or effective sample sizes.
real matrix rdistnn_tdnn_diagnostics(real matrix W, real matrix DX)
{
    real scalar j,total
    real colvector w
    real matrix out
    if(rows(W)!=rows(DX) | cols(W)!=cols(DX)) exit(error(198))
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
