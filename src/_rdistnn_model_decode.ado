version 17
// Version-1 numeric snapshot decoding. No score reconstruction or posting.
capture mata: mata drop rdistnn_model_decode()
capture mata: mata drop rdistnn_model_arms()
mata:
void rdistnn_model_decode(string scalar target_name, real scalar s,
    real scalar do_se, string scalar cmd, string scalar method,
    real scalar k, real matrix X, real matrix X0, real matrix signals,
    real matrix ranks, real matrix DX)
{
    string scalar distance
    real scalar n,m,q,j
    real matrix old_targets,geometry
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
        geometry=J(0,0,.)
        if(distance=="mahalanobis" | distance=="l1_scaled" | distance=="l2_scaled") {
            geometry=st_matrix("e(model_geometry)")
        }
        DX=rdistnn_distance_matrix(X,X0,distance,geometry)
        ranks=rdistnn_localization_ranks(DX,X)
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
}
void rdistnn_model_arms(string scalar cmd, real matrix A, real matrix E)
{
    string scalar arm,prop
    string rowvector vars
    real rowvector index
    real matrix data
    real scalar n,clip
    n=rows(st_matrix("e(model_X)"))
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
}
end
