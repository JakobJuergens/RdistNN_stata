version 17
// Version-3 numeric snapshot decoding. No score reconstruction or posting.
capture mata: mata drop rdistnn_model_invalid()
capture mata: mata drop rdistnn_model_validate()
capture mata: mata drop rdistnn_model_decode()
capture mata: mata drop rdistnn_model_arms()
mata:
// The numeric snapshot is authoritative only when its pieces agree. Validate
// cached geometry without refitting distances or modifying any e() value.
void rdistnn_model_invalid(string scalar reason)
{
    rdistnn_fail(498,"Saved model is invalid ("+rdistnn_message_value(reason)+"); refit with the current package.")
}

void rdistnn_model_validate()
{
    string scalar cmd,method,distance,field,se_method,expected_se
    string rowvector vars,xvars,required
    real scalar n,p,m,q,k,s,j,i,r,a,b,c,first,group,columns,clip,tol
    real matrix X,T,S,R,D,I,G,ids,groups,B,H,V,C
    real colvector ord,variance
    cmd=st_global("e(cmd)"); method=st_global("e(method)")
    distance=st_global("e(distance)")
    if(st_numscalar("e(model_version)")!=3 |
        st_global("e(tie_rule)")!="exact_distance_block_averaging_v1" |
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
    D=st_matrix("e(model_distances)"); groups=st_matrix("e(model_groups)"); I=st_matrix("e(model_inputs)")
    ids=st_matrix("e(model_row_ids)")
    n=rows(X); p=cols(X); m=rows(T); q=cols(S)
    vars=tokens(st_global("e(model_vars)")); xvars=tokens(st_global("e(xvars)"))
    if(n<1 | p<1 | m<1 | st_numscalar("e(N)")!=n |
        st_numscalar("e(p)")!=p | st_numscalar("e(kx)")!=p |
        st_numscalar("e(m)")!=m | cols(T)!=p | rows(S)!=n |
        rows(R)!=n | cols(R)!=m | rows(D)!=n | cols(D)!=m |
        rows(groups)!=n | cols(groups)!=m |
        rows(I)!=n | cols(I)!=cols(vars) | cols(vars)<1 | cols(xvars)!=p |
        rows(ids)!=n | cols(ids)!=1) rdistnn_model_invalid("snapshot dimensions")
    if(hasmissing(X) | hasmissing(T) | hasmissing(S) | hasmissing(I) |
        hasmissing(groups) | hasmissing(R) | hasmissing(D) | any(D:<0) | hasmissing(ids) |
        any(ids:<1) | any(ids:!=floor(ids))) rdistnn_model_invalid("nonfinite values, distances or observation IDs")
    // Original fitted observation IDs are strictly increasing, including when
    // if/in/markout selected a nonconsecutive subset. Do not bound them by n
    // or by the number of rows in the unrelated current dataset.
    if(n>1) if(any(ids[2..n]:<=ids[1..(n-1)])) rdistnn_model_invalid("observation IDs")
    if(q!=((cmd=="clate_dnn") ? 2 : ((cmd=="plqr_dnn") ? st_numscalar("e(kd)") : 1)) |
        q<1) rdistnn_model_invalid("signal components")
    // Cached reports must describe the same targets as reconstruction. Check
    // their shape and admissible uncertainty without recomputing estimates.
    C=st_matrix("e(targets)")
    if(rows(C)!=m | cols(C)!=p) rdistnn_model_invalid("report target dimensions")
    if(any(C:!=T)) rdistnn_model_invalid("report target identities")
    columns=(cmd=="plqr_dnn") ? m*q : m
    B=st_matrix("e(estimates)"); H=st_matrix("e(se_targets)")
    V=st_matrix("e(estimate_vcov)")
    if(rows(B)!=1 | cols(B)!=columns | rows(H)!=1 | cols(H)!=columns |
        rows(V)!=columns | cols(V)!=columns) rdistnn_model_invalid("canonical report dimensions")
    if(cmd!="clate_dnn" & hasmissing(B)) rdistnn_model_invalid("nonfinite estimates")
    if(any(H:<0)) rdistnn_model_invalid("negative standard errors")
    se_method=st_global("e(se_method)")
    expected_se=(cmd=="clate_dnn") ? "paired_rank_shift" :
        ((cmd=="np_tdnn") ? "signed_rank_shift" : "rank_shift")
    if(!any(se_method:==("none",expected_se)) |
        st_global("e(jackknife)")!=(se_method=="none" ? "none" : "fast")) {
        rdistnn_model_invalid("uncertainty method")
    }
    if(any(diagonal(V):<0) | any(V:!=V')) rdistnn_model_invalid("canonical covariance")
    if(any(H:!=sqrt(diagonal(V))')) rdistnn_model_invalid("standard error and covariance agreement")
    if(se_method=="none") {
        if(any(H:<.) | any(V:<.)) rdistnn_model_invalid("disabled uncertainty")
    }
    else if(cmd!="clate_dnn" & (hasmissing(H) | hasmissing(V))) {
        rdistnn_model_invalid("unavailable ordinary uncertainty")
    }
    G=st_matrix("e(diagnostics)")
    columns=(cmd=="np_tdnn") ? 8 : ((cmd=="np_dnn") ? 17 :
        ((cmd=="cate_dnn" | cmd=="clate_dnn") ? 25 : 11))
    if(rows(G)!=m | cols(G)!=columns) rdistnn_model_invalid("diagnostic dimensions")
    first=(cmd=="np_tdnn") ? 5 : 1
    if(hasmissing(G[,first]) | any(G[,first]:<=0)) {
        rdistnn_model_invalid("diagnostic localization concentration")
    }
    if((cmd=="cate_dnn" & st_global("e(score_method)")=="aipw") |
        (cmd=="clate_dnn" & st_global("e(score_method)")=="augmented_iv")) {
        clip=st_numscalar("e(clip)")
        if(missing(clip) | clip<=0 | clip>=.5 | 1-clip>=1) {
            rdistnn_model_invalid("propensity clipping control")
        }
        required=(cmd=="cate_dnn") ? ("t_col","e_col","m0_col","m1_col") :
            ("d_col","a_col","e_col","my0_col","my1_col","md1_col")
        for(j=1;j<=cols(required);j++) {
            field=st_global("e("+required[j]+")")
            if(field=="" | sum(vars:==field)!=1) rdistnn_model_invalid("retained nuisance input links")
        }
    }
    if(cmd=="clate_dnn") {
        tol=st_numscalar("e(denominator_tol)")
        if(missing(tol) | tol<0) rdistnn_model_invalid("denominator tolerance")
        C=st_matrix("e(status)")
        if(rows(C)!=m | cols(C)!=5) rdistnn_model_invalid("ratio status dimensions")
        if(hasmissing(C) | any(C:!=floor(C)) | any(C:<0) |
            any(C[,1]:>2) | any(C[,4]:>3) | any(C[,(2,3,5)]:>1)) {
            rdistnn_model_invalid("ratio status values")
        }
        if(any((B':>=.):!=(C[,1]:!=0))) rdistnn_model_invalid("ratio estimate availability")
        if(se_method=="none") {
            if(any(C[,4]:!=1) | any(C[,5]:!=0)) rdistnn_model_invalid("disabled ratio uncertainty status")
        }
        else {
            if(any(C[,4]:==1) |
                any(C[,5]:!=((C[,1]:==0):&(C[,4]:==0):&(diagonal(V):>0):&(diagonal(V):<.)))) {
                rdistnn_model_invalid("ratio uncertainty status")
            }
            variance=diagonal(V)
            if(any((C[,1]:==0):&(C[,4]:==0):&((variance:<=0):|(variance:>=.))) |
                any((C[,4]:==2):&((C[,1]:!=0):|(variance:!=0))) |
                any(((C[,4]:==3):|(C[,1]:!=0)):&(variance:<.))) {
                rdistnn_model_invalid("ratio covariance status and availability")
            }
            for(j=1;j<=m;j++) if(missing(variance[j])) {
                if(any(V[j,]:<.)) rdistnn_model_invalid("unavailable ratio covariance row")
            }
        }
        C=st_matrix("e(components)")
        if(rows(C)!=m | cols(C)!=2) rdistnn_model_invalid("paired component dimensions")
        C=st_matrix("e(component_vcov)")
        if(rows(C)!=2*m | cols(C)!=2*m) rdistnn_model_invalid("paired covariance dimensions")
    }
    k=st_numscalar("e(k)"); s=st_numscalar("e(s)")
    if(hasmissing((k,s)) | k<1 | k!=floor(k) | s<k | s>n | s!=floor(s)) {
        rdistnn_model_invalid("localization orders")
    }
    if(se_method!="none" & s>=n) rdistnn_model_invalid("delete-one uncertainty order")
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
    // Check group IDs in the same O(n*m) scan, without sorting or distance refits.
    if(any(R:<1) | any(R:>n) | any(R:!=floor(R))) rdistnn_model_invalid("rank permutation")
    for(j=1;j<=m;j++) {
        ord=J(n,1,0)
        for(i=1;i<=n;i++) {
            r=R[i,j]
            if(ord[r]!=0) rdistnn_model_invalid("rank permutation")
            ord[r]=i
        }
        group=1
        if(groups[ord[1],j]!=group) rdistnn_model_invalid("cached distance groups")
        for(i=2;i<=n;i++) {
            a=ord[i-1]; b=ord[i]
            if(D[a,j]>D[b,j]) rdistnn_model_invalid("rank-distance order")
            if(D[a,j]==D[b,j]) {
                if(ids[a]>ids[b]) rdistnn_model_invalid("distance tie observation order")
            }
            else group++
            if(groups[b,j]!=group) rdistnn_model_invalid("cached distance groups")
        }
    }
}

void rdistnn_model_decode(string scalar target_name, real scalar s,
    real scalar do_se, string scalar cmd, string scalar method,
    real scalar k, real matrix X, real matrix X0, real matrix signals,
    real matrix ranks, real matrix DX, real matrix groups)
{
    string scalar distance
    real scalar n,m,q,j,t
    real rowvector cached_columns
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
    // Exact retained target rows are the same geometry even when callers
    // repeat or reorder them. Validation above precedes every cached read.
    cached_columns=J(1,m,0)
    for(j=1;j<=m;j++) {
        for(t=1;t<=rows(old_targets);t++) {
            if(all(old_targets[t,]:==X0[j,])) {
                cached_columns[j]=t
                break
            }
        }
    }
    if(all(cached_columns:>0)) {
        ranks=st_matrix("e(model_ranks)")[,cached_columns]
        DX=st_matrix("e(model_distances)")[,cached_columns]
        groups=st_matrix("e(model_groups)")[,cached_columns]
    }
    if(rows(ranks)==0) {
        geometry=J(0,0,.)
        if(distance=="mahalanobis" | distance=="l1_scaled" | distance=="l2_scaled") {
            geometry=st_matrix("e(model_geometry)")
        }
        DX=rdistnn_distance_matrix(X,X0,distance,geometry)
        ranks=rdistnn_localization_ranks(DX)
        groups=rdistnn_distance_groups(ranks,DX,1)
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
