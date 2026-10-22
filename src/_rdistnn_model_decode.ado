version 17
// Version-2 numeric snapshot decoding. No score reconstruction or posting.
capture mata: mata drop rdistnn_model_invalid()
capture mata: mata drop rdistnn_model_validate()
capture mata: mata drop rdistnn_model_decode()
capture mata: mata drop rdistnn_model_arms()
mata:
// The numeric snapshot is authoritative only when its pieces agree. Validate
// cached geometry without refitting distances or modifying any e() value.
void rdistnn_model_invalid(string scalar reason)
{
    errprintf("Saved model is invalid (%s); refit with the current package.\n",reason)
    exit(error(498))
}

void rdistnn_model_validate()
{
    string scalar cmd,method,distance
    string rowvector vars,xvars
    real scalar n,p,m,q,k,s,j,i,r,a,b,c,first
    real matrix X,T,S,R,D,I,G,ids
    real colvector ord
    cmd=st_global("e(cmd)"); method=st_global("e(method)")
    distance=st_global("e(distance)")
    if(st_numscalar("e(model_version)")!=2 |
        !any(cmd:==("np_dnn","np_tdnn","cate_dnn","clate_dnn","plqr_dnn")) |
        !any(method:==("exact","asymp")) |
        !any(distance:==("l1","l2","l1_scaled","l2_scaled","mahalanobis"))) {
        rdistnn_model_invalid("schema, command, method or distance")
    }
    if(cmd=="cate_dnn" & !any(st_global("e(score_method)"):==("aipw","pseudo_outcome"))) {
        rdistnn_model_invalid("CATE score method")
    }
    if(cmd=="clate_dnn" & !any(st_global("e(score_method)"):==("augmented_iv","paired_signals"))) {
        rdistnn_model_invalid("CLATE score method")
    }
    X=st_matrix("e(model_X)"); T=st_matrix("e(model_targets)")
    S=st_matrix("e(model_signals)"); R=st_matrix("e(model_ranks)")
    D=st_matrix("e(model_distances)"); I=st_matrix("e(model_inputs)")
    ids=st_matrix("e(model_row_ids)")
    n=rows(X); p=cols(X); m=rows(T); q=cols(S)
    vars=tokens(st_global("e(model_vars)")); xvars=tokens(st_global("e(xvars)"))
    if(n<1 | p<1 | m<1 | st_numscalar("e(N)")!=n |
        st_numscalar("e(p)")!=p | st_numscalar("e(kx)")!=p |
        st_numscalar("e(m)")!=m | cols(T)!=p | rows(S)!=n |
        rows(R)!=n | cols(R)!=m | rows(D)!=n | cols(D)!=m |
        rows(I)!=n | cols(I)!=cols(vars) | cols(vars)<1 | cols(xvars)!=p |
        rows(ids)!=n | cols(ids)!=1) rdistnn_model_invalid("snapshot dimensions")
    if(hasmissing(X) | hasmissing(T) | hasmissing(S) | hasmissing(I) |
        hasmissing(R) | hasmissing(D) | any(D:<0) | hasmissing(ids) |
        any(ids:<1) | any(ids:!=floor(ids))) rdistnn_model_invalid("nonfinite values, distances or observation IDs")
    // Original fitted observation IDs are strictly increasing, including when
    // if/in/markout selected a nonconsecutive subset. Do not bound them by n
    // or by the number of rows in the unrelated current dataset.
    if(n>1) if(any(ids[2..n]:<=ids[1..(n-1)])) rdistnn_model_invalid("observation IDs")
    if(q!=((cmd=="clate_dnn") ? 2 : ((cmd=="plqr_dnn") ? st_numscalar("e(kd)") : 1)) |
        q<1) rdistnn_model_invalid("signal components")
    k=st_numscalar("e(k)"); s=st_numscalar("e(s)")
    if(hasmissing((k,s)) | k<1 | k!=floor(k) | s<k | s>n | s!=floor(s)) {
        rdistnn_model_invalid("localization orders")
    }
    if(cmd=="np_tdnn") {
        a=st_numscalar("e(s1)"); b=st_numscalar("e(s2)")
        c=st_numscalar("e(dim)"); first=st_numscalar("e(normalize)")
        if(hasmissing((a,b,c,first)) | a<1 | a!=floor(a) | b<=a |
            b!=floor(b) | b>n | b!=s | k!=1 | c<1 | c!=floor(c) |
            (first!=0 & first!=1)) rdistnn_model_invalid("two-scale controls")
    }
    if(distance=="l1_scaled" | distance=="l2_scaled") {
        G=st_matrix("e(model_geometry)")
        if(rows(G)!=1 | cols(G)!=p | hasmissing(G) | any(G:<=0)) {
            rdistnn_model_invalid("frozen scales")
        }
    }
    if(distance=="mahalanobis") {
        G=st_matrix("e(model_geometry)")
        if(rows(G)!=p+2 | cols(G)!=p | hasmissing(G)) rdistnn_model_invalid("frozen whitening geometry")
        if(G[1,1]<=0) rdistnn_model_invalid("frozen whitening scale")
    }
    // Invert each rank permutation once, then check adjacent ordered rows.
    // This is O(n*m*p) in tied coordinates, without sorting or distance refits.
    if(any(R:<1) | any(R:>n) | any(R:!=floor(R))) rdistnn_model_invalid("rank permutation")
    for(j=1;j<=m;j++) {
        ord=J(n,1,0)
        for(i=1;i<=n;i++) {
            r=R[i,j]
            if(ord[r]!=0) rdistnn_model_invalid("rank permutation")
            ord[r]=i
        }
        for(i=2;i<=n;i++) {
            a=ord[i-1]; b=ord[i]
            if(D[a,j]>D[b,j]) rdistnn_model_invalid("rank-distance order")
            if(D[a,j]!=D[b,j]) continue
            for(c=1;c<=p;c++) {
                if(X[a,c]<X[b,c]) break
                if(X[a,c]>X[b,c]) rdistnn_model_invalid("covariate tie order")
            }
        }
    }
}

void rdistnn_model_decode(string scalar target_name, real scalar s,
    real scalar do_se, string scalar cmd, string scalar method,
    real scalar k, real matrix X, real matrix X0, real matrix signals,
    real matrix ranks, real matrix DX)
{
    string scalar distance
    real scalar n,m,q
    real matrix old_targets,geometry
    rdistnn_model_validate()
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
