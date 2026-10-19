*! version 0.1.0  2026-10-18
version 17
// Signed two-scale rank tables and generic, frozen-rank deletion inference.
foreach fn in coefficients table weights loo_tables cov diagnostics {
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
        for(r=1;r<=n-s1+1;r++) a[r]=exp(lncomb(n-r,s1-1)-lncomb(n,s1))
        for(r=1;r<=n-s2+1;r++) b[r]=exp(lncomb(n-r,s2-1)-lncomb(n,s2))
        total=quadsum(a)
        if(total<=0 | missing(total)) exit(error(498))
        if(abs(total-1)>1e-12) a=a/total
        total=quadsum(b)
        if(total<=0 | missing(total)) exit(error(498))
        if(abs(total-1)>1e-12) b=b/total
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
        total=quadsum(w)
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

// Compatibility wrapper over the common arbitrary-table dispatcher.
real matrix rdistnn_tdnn_loo_tables(real matrix ranks, real matrix signal,
    real colvector wn, real colvector wnm1, real scalar center,
    | string scalar engine)
{
    if(rows(wn)!=rows(ranks) | cols(wn)!=1 | hasmissing(wn)) {
        errprintf("Full rank table must contain n finite entries.\n")
        exit(error(198))
    }
    if(args()<6) engine=""
    return(rdistnn_loo_from_table(ranks,signal,wnm1,center,engine))
}

real matrix rdistnn_tdnn_cov(real matrix ranks, real matrix signal,
    real scalar s1, real scalar s2, real scalar dim, string scalar method,
    real scalar normalize)
{
    real scalar n
    real matrix loo,V
    n=rows(ranks)
    loo=rdistnn_tdnn_loo_tables(ranks,signal,
        rdistnn_tdnn_table(n,s1,s2,dim,method,normalize),
        rdistnn_tdnn_table(n-1,s1,s2,dim,method,normalize),1)
    loo=(loo:-mean(loo))*sqrt((n-1)/n)
    V=cross(loo,loo)
    return((V+V')/2)
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
