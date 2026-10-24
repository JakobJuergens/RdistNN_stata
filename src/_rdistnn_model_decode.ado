version 17
// Version-4 numeric snapshot decoding. No score reconstruction or posting.
capture mata: mata drop rdistnn_model_invalid()
capture mata: mata drop rdistnn_model_validate()
capture mata: mata drop rdistnn_model_validate_grid()
capture mata: mata drop rdistnn_model_decode()
capture mata: mata drop rdistnn_model_arms()
capture mata: mata drop rdistnn_model_target_index()
mata:
// An empty public selection means every target. Keep that sequence numeric:
// Stata numlists have a lower size limit than supported result matrices.
real rowvector rdistnn_model_target_index(real rowvector targets, real scalar m)
{
    if(cols(targets)==0) return(1..m)
    return(targets)
}

// The numeric snapshot is authoritative only when its pieces agree. Validate
// cached geometry without refitting distances or modifying any e() value.
void rdistnn_model_invalid(string scalar reason)
{
    rdistnn_fail(498,"Saved model is invalid ("+rdistnn_message_value(reason)+"); refit with the current package.")
}

// Grid payloads are part of the current snapshot, not alternate fitted models.
// Check their identities and availability without reconstructing any weights.
void rdistnn_model_validate_grid(string scalar cmd, real scalar n,
    real scalar m, real scalar k, real scalar s, string scalar se_method)
{
    string scalar mode,source
    real scalar g,h,ng,default_id,reason,cid,nu,nuf
    real matrix C,E,H,D,S,Y,A,one
    real colvector seen,representative
    mode=st_global("e(tuning_mode)")
    if(cmd=="plqr_dnn" & mode=="") return
    if(!any(mode:==("default","grid"))) rdistnn_model_invalid("tuning mode")
    if(mode!="grid") return
    if(!any(cmd:==("np_dnn","cate_dnn","clate_dnn"))) rdistnn_model_invalid("candidate command")
    C=st_matrix("e(tuning_candidates)"); ng=rows(C)
    E=st_matrix("e(tuning_estimates)"); H=st_matrix("e(tuning_se)")
    D=st_matrix("e(tuning_diagnostics)")
    if(ng<1 | cols(C)!=8 | rows(E)!=ng | cols(E)!=m |
        rows(H)!=ng | cols(H)!=m | rows(D)!=ng | cols(D)!=5) {
        rdistnn_model_invalid("candidate table dimensions")
    }
    if(hasmissing(C[,1..6]) | any(C[,1]:!=(1::ng)) |
        any(C[,2]:<1) | any(C[,2]:!=floor(C[,2])) |
        any((C[,3]:!=0):&(C[,3]:!=1)) | sum(C[,3])!=1 |
        any(C[,7]:>.) | any((C[,7]:<.):&(C[,7]:<=0)) |
        any(C[,8]:>.) | any(C[,8]:<0) |
        any((C[,8]:<.):&(C[,8]:!=floor(C[,8])))) {
        rdistnn_model_invalid("candidate identities, flags or orders")
    }
    default_id=selectindex(C[,3])[1]
    if(st_numscalar("e(tuning_default_candidate)")!=default_id |
        C[default_id,8]!=s | C[default_id,4]!=1) rdistnn_model_invalid("candidate default anchor")
    source=st_global("e(tuning_source)")
    if(!any(source:==("sgrid","lambda_grid"))) rdistnn_model_invalid("candidate source")
    if(source=="sgrid") {
        if(any(C[,7]:!=.) | hasmissing(C[,8]) | any(C[,8]:<1) |
            st_numscalar("e(lambda)")!=.) rdistnn_model_invalid("direct candidate controls")
    }
    else if(missing(st_numscalar("e(lambda)")) | st_numscalar("e(lambda)")<=0) {
        rdistnn_model_invalid("candidate lambda control")
    }
    if(any(H:<0) | any(E[default_id,]:!=st_matrix("e(estimates)")) |
        any(H[default_id,]:!=st_matrix("e(se_targets)"))) rdistnn_model_invalid("candidate default estimates or uncertainty")
    if(any(st_matrixcolstripe("e(tuning_estimates)"):!=st_matrixcolstripe("e(estimates)")) |
        any(st_matrixcolstripe("e(tuning_se)"):!=st_matrixcolstripe("e(estimates)"))) {
        rdistnn_model_invalid("candidate target identities")
    }
    if(cmd=="clate_dnn") {
        Y=st_matrix("e(tuning_reduced_form)"); A=st_matrix("e(tuning_first_stage)")
        S=st_matrix("e(tuning_target_status)")
        if(rows(Y)!=ng | cols(Y)!=m | rows(A)!=ng | cols(A)!=m |
            rows(S)!=ng*m | cols(S)!=7) rdistnn_model_invalid("paired candidate dimensions")
        if(any(Y[default_id,]:!=st_matrix("e(components)")[,1]') |
            any(A[default_id,]:!=st_matrix("e(components)")[,2]') |
            any(S[(default_id-1)*m:+(1::m),3..7]:!=st_matrix("e(status)"))) {
            rdistnn_model_invalid("paired candidate default anchor")
        }
    }
    seen=representative=J(0,1,.); nu=nuf=0
    for(g=1;g<=ng;g++) {
        reason=0
        if(k>C[g,8]) reason=1
        else if(se_method!="none" & C[g,8]>=n) reason=2
        else if(se_method=="none" & C[g,8]>n) reason=3
        if(C[g,4]!=(reason==0) | C[g,5]!=reason | C[g,6]!=C[g,4]) {
            rdistnn_model_invalid("candidate feasibility or status")
        }
        cid=C[g,2]
        if(cid>nu+1) rdistnn_model_invalid("candidate computed-order identities")
        if(cid==nu+1) {
            if(any(seen:==C[g,8])) rdistnn_model_invalid("candidate duplicate-order identity")
            seen=seen\C[g,8]; representative=representative\g
            nu++; nuf=nuf+C[g,4]
        }
        else {
            h=representative[cid]
            if(C[g,8]!=seen[cid] | any(E[g,]:!=E[h,]) |
                any(H[g,]:!=H[h,]) | any(D[g,]:!=D[h,])) rdistnn_model_invalid("candidate alias payload")
            if(cmd=="clate_dnn") {
                if(any(Y[g,]:!=Y[h,]) | any(A[g,]:!=A[h,]) |
                    any(S[(g-1)*m:+(1::m),2..7]:!=S[(h-1)*m:+(1::m),2..7])) {
                    rdistnn_model_invalid("paired candidate alias payload")
                }
            }
        }
        if(D[g,3]!=min((5,n))) rdistnn_model_invalid("candidate diagnostic support count")
        if(C[g,4]==0) {
            if(any(E[g,]:<.) | any(H[g,]:<.) | any(D[g,(1,2,4,5)]:<.)) {
                rdistnn_model_invalid("infeasible candidate availability")
            }
            if(cmd=="clate_dnn") {
                if(any(Y[g,]:<.) | any(A[g,]:<.) | any(S[(g-1)*m:+(1::m),]:<.)) {
                    rdistnn_model_invalid("infeasible paired candidate availability")
                }
            }
            continue
        }
        if(hasmissing(D[g,]) | D[g,1]<=0 | D[g,2]<0 |
            D[g,4]<0 | D[g,4]>1+8*epsilon(1) | (D[g,5]!=0 & D[g,5]!=1)) {
            rdistnn_model_invalid("candidate diagnostic values")
        }
        if(cmd!="clate_dnn" & hasmissing(E[g,])) rdistnn_model_invalid("candidate estimate availability")
        if(se_method=="none") {
            if(any(H[g,]:<.)) rdistnn_model_invalid("disabled candidate uncertainty")
        }
        else if(cmd!="clate_dnn" & hasmissing(H[g,])) rdistnn_model_invalid("candidate uncertainty availability")
        if(cmd=="clate_dnn") {
            one=S[(g-1)*m:+(1::m),]
            if(hasmissing(one) |
                any(one[,1]:!=g) | any(one[,2]:!=(1::m)) |
                any(one[,3..7]:!=floor(one[,3..7])) | any(one[,3..7]:<0) |
                any(one[,3]:>2) | any(one[,6]:>3) | any(one[,(4,5,7)]:>1) |
                any((E[g,]':>=.):!=(one[,3]:!=0)) |
                any(((Y[g,]':>=.):|(A[g,]':>=.)):&(one[,3]:!=2))) rdistnn_model_invalid("paired candidate status identities or availability")
            if(se_method=="none") {
                if(any(one[,6]:!=1) | any(one[,7]:!=0)) rdistnn_model_invalid("disabled paired candidate status")
            }
            else if(any(one[,6]:==1) |
                any((one[,3]:==0):&(one[,6]:==0):&((H[g,]':<=0):|(H[g,]':>=.))) |
                any((one[,6]:==2):&((one[,3]:!=0):|(H[g,]':!=0))) |
                any(((one[,6]:==3):|(one[,3]:!=0)):&(H[g,]':<.)) |
                any(one[,7]:!=((one[,3]:==0):&(one[,6]:==0)))) {
                rdistnn_model_invalid("paired candidate uncertainty status")
            }
        }
    }
    if(st_numscalar("e(tuning_n_candidates)")!=ng |
        st_numscalar("e(tuning_n_feasible)")!=sum(C[,4]) |
        st_numscalar("e(tuning_n_skipped)")!=ng-sum(C[,4]) |
        st_numscalar("e(tuning_n_failed)")!=0 |
        st_numscalar("e(tuning_n_unique)")!=nuf) rdistnn_model_invalid("candidate counts")
}

void rdistnn_model_validate()
{
    string scalar cmd,method,distance,field,se_method,expected_se
    string rowvector vars,xvars,required,signature_vars
    real scalar n,p,m,q,k,s,j,i,r,a,b,c,first,group,columns,clip,tol
    real matrix X,T,S,R,D,I,G,ids,groups,B,H,V,C,notices
    real colvector ord,variance,ties,weak
    cmd=st_global("e(cmd)"); method=st_global("e(method)")
    distance=st_global("e(distance)")
    if(st_numscalar("e(model_version)")!=4 |
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
    signature_vars=tokens(st_global("e(signature_vars)"))
    if(cols(signature_vars)!=cols(vars) | cols(signature_vars)<1 |
        st_global("e(ordered_signature)")=="") rdistnn_model_invalid("retained input signature metadata")
    if(any(signature_vars:!=vars)) rdistnn_model_invalid("retained signature variable identities")
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
    notices=st_matrix("e(target_notices)")
    if(rows(notices)!=m | cols(notices)!=2) rdistnn_model_invalid("target notice dimensions")
    if(hasmissing(notices) | any((notices:!=0):&(notices:!=1))) rdistnn_model_invalid("target notice flags")
    ties=J(m,1,0)
    first=(cmd=="np_tdnn") ? 5 : 1
    if(hasmissing(G[,first]) | any(G[,first]:<=0)) {
        rdistnn_model_invalid("diagnostic localization concentration")
    }
    weak=(n>=2):&(G[,first]:>=.9*n)
    if(any(notices[,2]:!=weak)) rdistnn_model_invalid("weak-localization notice facts")
    if(st_numscalar("e(n_tied_targets)")!=sum(notices[,1]) |
        st_numscalar("e(n_weak_targets)")!=sum(weak)) rdistnn_model_invalid("target notice counts")
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
    rdistnn_model_validate_grid(cmd,n,m,k,s,se_method)
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
        ties[j]=(group<n)
    }
    if(any(notices[,1]:!=ties)) rdistnn_model_invalid("distance-tie notice facts")
}

void rdistnn_model_decode(string scalar target_name, real scalar s,
    real scalar do_se, string scalar cmd, string scalar method,
    real scalar k, real matrix X, real matrix X0, real matrix signals,
    real matrix ranks, real matrix DX, real matrix groups, | real scalar validated)
{
    string scalar distance
    real scalar n,m,q,j,t
    real rowvector cached_columns
    real matrix old_targets,geometry
    if(args()<13) validated=0
    if(!validated) rdistnn_model_validate()
    cmd=st_global("e(cmd)"); method=st_global("e(method)")
    distance=st_global("e(distance)"); k=st_numscalar("e(k)")
    X=st_matrix("e(model_X)"); signals=st_matrix("e(model_signals)")
    X0=st_matrix(target_name); n=rows(X); m=rows(X0); q=cols(signals)
    if(n<1 | m<1 | cols(X0)!=cols(X) | rows(signals)!=n |
        hasmissing(X) | hasmissing(X0) | hasmissing(signals) |
        s<k | s>n | s!=floor(s) | (do_se & s>=n)) {
        rdistnn_fail(198,"Saved model and targets require finite aligned matrices and k <= s <= n (s < n for covariance).")
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
        rdistnn_fail(498,"Saved geometry is incompatible with the fitted observations and targets; refit the model.")
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
