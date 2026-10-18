// Re-evaluate a saved model without reading the caller's dataset or retraining.
capture mata: mata drop rdistnn_model_diagnostics()
capture mata: mata drop rdistnn_model_evaluate()
capture mata: mata drop rdistnn_model_weights()
capture mata: mata drop rdistnn_model_plot_data()
mata:
real matrix rdistnn_model_plot_data(string scalar type, real rowvector targets)
{
    real matrix data,out,signals,one
    real colvector w,cum,ord,score,contribution
    real scalar j,t,n,paired
    data=rdistnn_model_weights(targets); n=rows(st_matrix("e(model_X)"))
    signals=st_matrix("e(model_signals)"); paired=st_global("e(cmd)")=="clate_dnn"
    out=J(0, type=="kernel" ? 5 : (paired ? 9 : 7),.)
    for(j=1;j<=cols(targets);j++) {
        t=targets[j]; one=data[(j-1)*n:+(1::n),]; w=one[,5]
        if(type=="kernel") {
            ord=order((one[,3],(1::n)),(1,2)); cum=J(n,1,.)
            cum[ord]=runningsum(w[ord]); out=out\(one[,2..5],cum)
        }
        else {
            if(missing(st_matrix("e(estimates)")[1,t])) {
                errprintf("Influence is unavailable for undefined target %g.\n",t)
                exit(error(498))
            }
            score=signals[,1]
            if(paired) score=score-st_matrix("e(estimates)")[1,t]*signals[,2]
            contribution=w:*(score:-(paired ? 0 : st_matrix("e(estimates)")[1,t]))
            one=(one[,2],one[,1],one[,3..5],score,contribution)
            if(paired) one=(one,signals)
            out=out\one
        }
    }
    return(out)
}

real matrix rdistnn_model_weights(real rowvector targets)
{
    real matrix R,D,W,out
    real scalar j,t,n
    R=st_matrix("e(model_ranks)"); D=st_matrix("e(model_distances)")
    W=rdistnn_weights_from_rmat(R,st_numscalar("e(k)"),st_numscalar("e(s)"),st_global("e(method)"))
    n=rows(R); out=J(0,5,.)
    for(j=1;j<=cols(targets);j++) {
        t=targets[j]
        out=out\(st_matrix("e(model_row_ids)"),J(n,1,t),R[,t],D[,t],W[,t])
    }
    return(out)
}

real matrix rdistnn_model_diagnostics(string scalar cmd, real matrix X,
    real matrix X0, real matrix signals, real matrix W, real matrix DX,
    real scalar k, real scalar s, real rowvector b, real matrix A, real matrix E)
{
    real scalar n,m,j,total,ybar
    real matrix out,pair
    real colvector w,d,ord,sorted
    n=rows(X); m=rows(X0)
    out=J(m,25,.)
    if(cmd=="np_dnn") out=J(m,17,.)
    if(cmd=="plqr_dnn") out=J(m,11,.)
    for(j=1;j<=m;j++) {
        w=W[,j]; d=DX[,j]; total=quadsum(w)
        ord=order((w,(1::n)),(-1,2)); sorted=sort(d,1)
        out[j,1..11]=(1/quadsum(w:^2),max(w),w[ord[1]],
            quadsum(w[ord[1..min((5,n))]]),quadsum(w[ord[1..min((10,n))]]),
            sum(abs(w):>sqrt(epsilon(1))),1,sorted[min((k,n))],
            sorted[min((s,n))],cross(d,w)/total,rdistnn_weighted_quantile(d,w,.5))
        if(cmd=="np_dnn") {
            ybar=cross(signals,w)/total
            out[j,12..17]=(1,ybar,sqrt(cross((signals:-ybar):^2,w)/total),
                rdistnn_weighted_quantile(signals,w,.1),
                rdistnn_weighted_quantile(signals,w,.5),
                rdistnn_weighted_quantile(signals,w,.9))
        }
        else if(cmd!="plqr_dnn") {
            out[j,20]=rdistnn_weighted_quantile(sqrt(rowsum((X:-X0[j,]):^2)),w,.9)
        }
    }
    if(cmd=="cate_dnn" | cmd=="clate_dnn") {
        pair=signals
        if(cmd=="cate_dnn") pair=(signals,J(n,1,1))
        out=rdistnn_clate_diagnostics(out,W,pair,b,A,E)
    }
    return(out)
}

// Output names: estimates, covariance, diagnostics, components, component V,
// status, effective sample sizes, ranks, distances, weights and q-by-m fits.
void rdistnn_model_evaluate(string scalar target_name, real scalar s,
    real scalar do_se, string rowvector names)
{
    string scalar cmd,method,distance,arm,prop
    string rowvector vars
    real scalar n,m,q,k,tol,j,clip
    real rowvector index
    real matrix X,X0,signals,ranks,DX,W,diag,components,CV,status,V,fit,A,E,data
    real matrix old_targets,loo,centered
    real rowvector b,neff
    cmd=st_global("e(cmd)"); method=st_global("e(method)")
    distance=st_global("e(distance)"); k=st_numscalar("e(k)")
    X=st_matrix("e(model_X)"); signals=st_matrix("e(model_signals)")
    X0=st_matrix(target_name); n=rows(X); m=rows(X0); q=cols(signals)
    if(n<1 | m<1 | cols(X0)!=cols(X) | rows(signals)!=n |
        hasmissing(X) | hasmissing(X0) | hasmissing(signals) |
        s<k | s>n | s!=floor(s) | (do_se & s>=n)) {
        errprintf("Saved model and targets require finite aligned matrices and k <= s <= n (s < n for covariance).\n")
        exit(error(198))
    }
    ranks=DX=J(0,0,.)
    old_targets=st_matrix("e(model_targets)")
    if(rows(old_targets)==m & cols(old_targets)==cols(X0)) {
        if(all(old_targets:==X0)) {
            ranks=st_matrix("e(model_ranks)")
            DX=st_matrix("e(model_distances)")
        }
    }
    if(rows(ranks)==0) {
        DX=rdistnn_distance_matrix(X,X0,distance)
        ranks=rdistnn_ranks_from_dmat(DX)
    }
    if(rows(ranks)!=n | cols(ranks)!=m | rows(DX)!=n | cols(DX)!=m | hasmissing(DX)) {
        errprintf("Saved geometry is incompatible with the fitted observations and targets.\n")
        exit(error(498))
    }
    for(j=1;j<=m;j++) {
        if(any(sort(ranks[,j],1):!=(1::n))) {
            errprintf("Saved ranks must be permutations of 1..n.\n")
            exit(error(498))
        }
    }
    W=rdistnn_weights_from_rmat(ranks,k,s,method)
    fit=signals'*W; b=vec(fit)'; neff=1:/colsum(W:^2)
    V=J(m*q,m*q,.)
    components=CV=status=J(1,1,.)
    if(cmd=="clate_dnn") {
        tol=st_numscalar("e(denominator_tol)")
        rdistnn_clate_fit(ranks,signals,k,s,method,do_se,tol,components,b,V,CV,status)
        fit=b
    }
    else if(do_se) {
        if(q==1) V=rdistnn_jk_cov_multi_from_ranks(ranks,signals,k,s,method)
        else {
            loo=rdistnn_loo_from_ranks(ranks,signals,k,s,method,1)
            centered=(loo:-mean(loo))*sqrt((n-1)/n)
            V=cross(centered,centered)
        }
    }
    A=E=J(0,1,.)
    arm=st_global("e(t_col)")
    if(cmd=="clate_dnn") arm=st_global("e(a_col)")
    prop=st_global("e(e_col)")
    data=st_matrix("e(model_inputs)"); vars=tokens(st_global("e(model_vars)"))
    if(arm!="") {
        index=selectindex(vars:==arm)
        if(length(index)==1) A=data[,index]
    }
    if(prop!="") {
        index=selectindex(vars:==prop)
        if(length(index)==1) {
            E=data[,index]; clip=st_numscalar("e(clip)")
            if(clip<.) E=rowmax((J(n,1,clip),rowmin((E,J(n,1,1-clip)))))
        }
    }
    diag=rdistnn_model_diagnostics(cmd,X,X0,signals,W,DX,k,s,b,A,E)
    st_matrix(names[1],b); st_matrix(names[2],V); st_matrix(names[3],diag)
    st_matrix(names[4],components); st_matrix(names[5],CV); st_matrix(names[6],status)
    st_matrix(names[7],neff); st_matrix(names[8],ranks); st_matrix(names[9],DX)
    st_matrix(names[10],W); st_matrix(names[11],fit)
}
end
