*! version 0.1.5  2026-10-25
version 17
// Re-evaluate a saved model without reading the caller's dataset or retraining.
// Mata resolves all referenced adapters when compiling the model functions.
findfile "_rdistnn_tdnn_mata.ado"
quietly run "`r(fn)'"
findfile "_rdistnn_model_decode.ado"
quietly run "`r(fn)'"
capture mata: mata drop rdistnn_model_evaluate()
capture mata: mata drop rdistnn_model_weights()
capture mata: mata drop rdistnn_model_plot_data()
mata:
real matrix rdistnn_model_plot_data(string scalar type, real rowvector targets,
    | real scalar validated)
{
    real matrix data,out,signals,one
    real colvector w,cum,ord,score,contribution,active,right
    real scalar j,t,n,paired
    if(args()<3) validated=0
    data=rdistnn_model_weights(targets,validated); n=rows(st_matrix("e(model_X)"))
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
                rdistnn_fail(498,"Influence is unavailable for undefined target "+strofreal(t)+". Inspect e(status) and e(components).")
            }
            score=signals[,1]
            if(paired) score=score-st_matrix("e(estimates)")[1,t]*signals[,2]
            right=J(n,1,st_matrix("e(estimates)")[1,t])
            if(paired) {
                right=J(n,1,0); active=selectindex(w:>0)
                if(length(active)) right[active]=st_matrix("e(estimates)")[1,t]*signals[active,2]
            }
            rdistnn_score_share(w,signals[,1],right,contribution)
            one=(one[,2],one[,1],one[,3..5],score,contribution)
            if(paired) one=(one,signals)
            out=out\one
        }
    }
    return(out)
}

real matrix rdistnn_model_weights(real rowvector targets, | real scalar validated)
{
    real matrix R,D,G,W,out,ids
    real scalar j,t,n
    if(args()<2) validated=0
    if(!validated) rdistnn_model_validate()
    R=st_matrix("e(model_ranks)")[,targets]
    D=st_matrix("e(model_distances)")[,targets]
    G=st_matrix("e(model_groups)")[,targets]
    ids=st_matrix("e(model_row_ids)")
    if(st_global("e(cmd)")=="np_tdnn") {
        W=rdistnn_tdnn_weights(R,st_numscalar("e(s1)"),st_numscalar("e(s2)"),
            st_numscalar("e(dim)"),st_global("e(method)"),st_numscalar("e(normalize)"),D,G)
    }
    else W=rdistnn_weights_from_rmat(R,st_numscalar("e(k)"),st_numscalar("e(s)"),st_global("e(method)"),D,G)
    n=rows(R); out=J(n*cols(targets),5,.)
    for(j=1;j<=cols(targets);j++) {
        t=targets[j]
        out[(j-1)*n:+(1::n),]=(ids,J(n,1,t),R[,j],D[,j],W[,j])
    }
    return(out)
}

// Output names: estimates, covariance, diagnostics, components, component V,
// status, effective sample sizes, ranks, distances, weights and q-by-m fits.
void rdistnn_model_evaluate(string scalar target_name, real scalar s,
    real scalar do_se, string rowvector names, | real scalar validated)
{
    string scalar cmd,method
    real scalar n,m,q,k,tol
    real matrix X,X0,signals,ranks,DX,G,W,diag,components,CV,status,V,fit,A,E,loo
    real rowvector b,neff
    if(args()<5) validated=0
    rdistnn_model_decode(target_name,s,do_se,cmd,method,k,X,X0,signals,ranks,DX,G,validated)
    n=rows(X); m=rows(X0); q=cols(signals)
    if(cmd=="np_tdnn") {
        W=rdistnn_tdnn_weights(ranks,st_numscalar("e(s1)"),st_numscalar("e(s2)"),
            st_numscalar("e(dim)"),method,st_numscalar("e(normalize)"),DX,G)
    }
    else W=rdistnn_weights_from_rmat(ranks,k,s,method,DX,G)
    fit=rdistnn_weighted_signals(signals,W); b=vec(fit)'; neff=1:/colsum(W:^2)
    V=J(m*q,m*q,.)
    components=CV=status=J(1,1,.)
    if(cmd=="clate_dnn") {
        tol=st_numscalar("e(denominator_tol)")
        rdistnn_clate_fit(ranks,signals,k,s,method,do_se,tol,components,b,V,CV,status,DX,G)
        fit=b
    }
    else if(do_se & cmd=="np_tdnn") {
        V=rdistnn_tdnn_cov(ranks,signals,st_numscalar("e(s1)"),st_numscalar("e(s2)"),
            st_numscalar("e(dim)"),method,st_numscalar("e(normalize)"),DX,G)
    }
    else if(do_se) {
        if(q==1) V=rdistnn_jk_cov_multi_from_ranks(ranks,signals,k,s,method,DX,G)
        else {
            loo=rdistnn_loo_from_ranks(ranks,signals,k,s,method,1,
                rdistnn_weights_from_ranks((1::(n-1)),k,s,method),DX,G)
            V=rdistnn_covariance(loo,(n-1)/n,1)
        }
    }
    rdistnn_model_arms(cmd,A,E)
    if(cmd=="np_tdnn") diag=rdistnn_tdnn_diagnostics(W,DX)
    else diag=rdistnn_diagnostics(cmd,X,X0,signals,W,DX,k,s,b,A,E)
    st_matrix(names[1],b); st_matrix(names[2],V); st_matrix(names[3],diag)
    st_matrix(names[4],components); st_matrix(names[5],CV); st_matrix(names[6],status)
    st_matrix(names[7],neff); st_matrix(names[8],ranks); st_matrix(names[9],DX)
    st_matrix(names[12],G)
    st_matrix(names[10],W); st_matrix(names[11],fit)
}
end

capture mata: mata drop rdistnn_model_tdnn_plot()
mata:
real matrix rdistnn_model_tdnn_plot(real rowvector targets)
{
    real matrix data
    data=rdistnn_model_weights(targets,1)
    data=data[order(data,(2,3)),]
    return(data[,(3,5,2)])
}
end
