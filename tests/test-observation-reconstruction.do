// Adversarial postestimation audit. Run from the repository root after adding
// src/ to adopath. The parent audit serializes licensed Stata execution.
// Presence is checked through r(matrices), never assignment _rc: an absent
// r(matrix) can be assigned successfully as a missing 1-by-1 matrix in Stata.
version 17
do "tests/test_framework.do"
capture mkdir "build"
capture mkdir "build/premerge"

do "tests/helpers/postestimation-state.do"

// Keep nested Mata control flow outside a Stata if-block that can be skipped:
// Stata's skipped-block brace parser otherwise treats Mata braces as its own.
capture mata: mata drop pmpe_influence_oracle()
mata:
void pmpe_influence_oracle()
{
    real scalar j,t,n,error
    real matrix influence,fit,iblock,kblock
    real colvector ids,y,d,score,centered
    influence=st_matrix("PMPE_I")
    ids=st_data(.,"source_id","original_sample")
    y=st_data(.,"y","original_sample")
    d=st_data(.,"denominator","original_sample")
    fit=st_matrix("PMPE_E")
    n=rows(ids)
    error=0
    for(j=1;j<=2;j++) {
        t=3-j
        iblock=influence[((j-1)*n+1)..(j*n),]
        kblock=st_matrix("PMPE_K")[((j-1)*n+1)..(j*n),]
        score=y
        centered=y:-fit[1,t]
        if(st_global("e(cmd)")=="clate_dnn") {
            score=centered=y:-fit[1,t]:*d
            error=max((error,max(abs(iblock[,8..9]-(y,d)))))
        }
        error=max((error,max(abs(iblock[,2]-ids)),max(abs(iblock[,3..5]-kblock[,2..4])),max(abs(iblock[,6]-score)),max(abs(iblock[,7]-iblock[,5]:*centered)),abs(sum(iblock[,7]))))
    }
    st_numscalar("PMPE_influence_ok",error<1e-10)
}
end

test_section "premerge postestimation: original sample geometry and frozen signals"
foreach family in cate np clate {
    rdnn_dgp_basic, n(48) seed(733)
    gen long source_id = _n
    gen double denominator = .7+.03*x1
    replace y = . in 13
    replace x1 = . in 17
    matrix PMPE_X0 = (-.4,0 \ .8,.2)
    local opts ""
    if "`family'"=="cate" local opts "pseudooutcome"
    if "`family'"=="clate" local opts "pseudooutcome denominator(denominator)"
    quietly `family'_dnn y if source_id>=5 & source_id<=43 in 3/46, ///
        xcols(x1 x2) x0matrix(PMPE_X0) s(6) `opts'
    gen byte original_sample = e(sample)
    local n = e(N)
    matrix PMPE_E = e(estimates)
    tempfile restricted
    quietly save "`restricted'", replace
    _pmpe_snapshot
    quietly `family'_dnn_plot, type(weights) view(kernel) targets(2 1 2) ///
        name(pmpe_kernel) replace nodraw
    local matrices : r(matrices)
    local present : list posof "plot_data" in matrices
    local nr = r(n_rows)
    matrix PMPE_K = r(plot_data)
    test_assert "[premerge `family'] if/in and complete-case sample determine kernel size" `present'>0 & `n'==37 & `nr'==2*`n' & rowsof(PMPE_K)==2*`n' & colsof(PMPE_K)==5 & PMPE_K[1,1]==2 & PMPE_K[`=`n'+1',1]==1
    // Independent geometry and k=1 combinatorics: C(n-r,s-1)/C(n,s).
    // The oracle deliberately does not call any rdistnn rank/weight helper.
    mata {
        PMPE_X = st_data(.,"x1 x2","original_sample")
        PMPE_out = st_matrix("PMPE_K")
        PMPE_targets = st_matrix("PMPE_X0")
        PMPE_werror = PMPE_derror = PMPE_rerror = PMPE_cerror = 0
        for (PMPE_j=1; PMPE_j<=2; PMPE_j++) {
            PMPE_t = (3-PMPE_j)
            PMPE_dist = sqrt(rowsum((PMPE_X:-PMPE_targets[PMPE_t,]):^2))
            PMPE_order = order((PMPE_dist,(1::rows(PMPE_X))),(1,2))
            PMPE_ranks = J(rows(PMPE_X),1,.)
            PMPE_ranks[PMPE_order] = (1::rows(PMPE_X))
            PMPE_weights = J(rows(PMPE_X),1,0)
            for (PMPE_i=1; PMPE_i<=rows(PMPE_X); PMPE_i++) {
                PMPE_r = PMPE_ranks[PMPE_i]
                if (PMPE_r<=rows(PMPE_X)-5) PMPE_weights[PMPE_i] = exp(lnfactorial(rows(PMPE_X)-PMPE_r)-lnfactorial(5)-lnfactorial(rows(PMPE_X)-PMPE_r-5)-lnfactorial(rows(PMPE_X))+lnfactorial(6)+lnfactorial(rows(PMPE_X)-6))
            }
            PMPE_cum = J(rows(PMPE_X),1,.)
            PMPE_cum[PMPE_order] = runningsum(PMPE_weights[PMPE_order])
            PMPE_block = PMPE_out[((PMPE_j-1)*rows(PMPE_X)+1)..(PMPE_j*rows(PMPE_X)),]
            PMPE_derror = max((PMPE_derror,max(abs(PMPE_block[,3]-PMPE_dist))))
            PMPE_rerror = max((PMPE_rerror,max(abs(PMPE_block[,2]-PMPE_ranks))))
            PMPE_werror = max((PMPE_werror,max(abs(PMPE_block[,4]-PMPE_weights))))
            PMPE_cerror = max((PMPE_cerror,max(abs(PMPE_block[,5]-PMPE_cum))))
        }
        st_numscalar("PMPE_geometry_ok", PMPE_derror<1e-12 & PMPE_rerror==0 & PMPE_werror<1e-12 & PMPE_cerror<1e-12)
    }
    test_assert "[premerge `family'] kernel distances/ranks/weights/cumulative shares match independent metric distance and observation-ID oracle" scalar(PMPE_geometry_ok)==1
    if "`family'"!="np" {
        quietly `family'_dnn_plot, type(influence) targets(2 1 2) nlabel(0) ///
            name(pmpe_influence) replace nodraw
        local nr = r(n_rows)
        local matrices : r(matrices)
        local present : list posof "plot_data" in matrices
        matrix PMPE_I = r(plot_data)
        local icols = cond("`family'"=="clate",9,7)
        test_assert "[premerge `family'] influence has exact selected sample rows and family schema" `present'>0 & `nr'==2*`n' & rowsof(PMPE_I)==2*`n' & colsof(PMPE_I)==`icols'
        mata: pmpe_influence_oracle()
        test_assert "[premerge `family'] influence preserves original dataset IDs and exact component/score/weight identities" scalar(PMPE_influence_ok)==1
    }
    _pmpe_unchanged, data("`restricted'")
    local unchanged = r(ok)
    test_assert "[premerge `family'] observation views preserve all original variables and e()" `unchanged'==1

    // Edits outside the actual estimation sample must not invalidate a fit.
    replace y = y+123 in 1
    replace x1 = x1+999 in 2
    quietly `family'_dnn_plot, type(weights) view(kernel) targets(2 1) name(pmpe_outside) replace nodraw
    matrix PMPE_OUTSIDE = r(plot_data)
    mata: st_numscalar("PMPE_error",all(st_matrix("PMPE_OUTSIDE"):==st_matrix("PMPE_K")))
    test_assert "[premerge `family'] excluded-row edits leave sample geometry usable and unchanged" scalar(PMPE_error)==1
    // A missingness repair does not retroactively add an observation to e(sample).
    replace y = 1 in 13
    capture noisily `family'_dnn_plot, type(weights) view(kernel) name(pmpe_repair) replace nodraw
    local repair_rc = _rc
    test_assert "[premerge `family'] repaired formerly marked-out row does not change stored estimation sample" `repair_rc'==0
    gsort -source_id
    quietly `family'_dnn_plot, type(estimates) name(pmpe_stored) replace nodraw
    local matrices : r(matrices)
    local present : list posof "plot_data" in matrices
    matrix PMPE_STORED = r(plot_data)
    test_assert "[premerge `family'] stored estimate view remains independent of changed dataset order" `present'>0 & rowsof(PMPE_STORED)==2 & colsof(PMPE_STORED)==10
}

test_section "premerge postestimation: supplied nuisance signals remain frozen"
foreach family in cate clate {
    rdnn_dgp_basic, n(40) seed(380)
    gen byte a = D
    gen byte uptake = a
    gen double md1 = .8
    gen double md0 = 0
    gen double phi_y = m1-m0 + a*(y-m1)/e - (1-a)*(y-m0)/(1-e)
    gen double phi_d = md1-md0 + a*(uptake-md1)/e - (1-a)*(uptake-md0)/(1-e)
    if "`family'"=="cate" quietly cate_dnn y, xcols(x1 x2) x0(0 0) ///
        tcol(a) ecol(e) m0col(m0) m1col(m1) s(5)
    else quietly clate_dnn y, xcols(x1 x2) x0(0 0) dcol(uptake) acol(a) ///
        ecol(e) my0col(m0) my1col(m1) md0col(md0) md1col(md1) s(5)
    matrix PMPE_E = e(estimates)
    quietly `family'_dnn_plot, type(influence) name(pmpe_nuisance) replace nodraw
    matrix PMPE_NI = r(plot_data)
    mata {
        PMPE_ni=st_matrix("PMPE_NI")
        PMPE_py=st_data(.,"phi_y")
        PMPE_pd=st_data(.,"phi_d")
        PMPE_nierror=max(abs(PMPE_ni[,6]-PMPE_py))
        if (st_global("e(cmd)")=="clate_dnn") PMPE_nierror=max((max(abs(PMPE_ni[,8..9]-(PMPE_py,PMPE_pd))), max(abs(PMPE_ni[,6]-(PMPE_py:-st_matrix("PMPE_E")[1,1]:*PMPE_pd)))))
        st_numscalar("PMPE_error",PMPE_nierror)
    }
    test_assert "[premerge `family'] nuisance-route influence reconstructs the independently supplied signals" scalar(PMPE_error)<1e-12
    // Derived variables never supplied to the fit may change freely.
    replace phi_y = 1e6
    capture noisily `family'_dnn_plot, type(influence) name(pmpe_unused) replace nodraw
    local derived_rc = _rc
    test_assert "[premerge `family'] unrelated derived signal is not part of supplied-nuisance signature" `derived_rc'==0
    // Fitted numeric observations and frozen signals remain authoritative.
    replace m0 = m0+.001 in 1
    quietly `family'_dnn_plot, type(influence) name(pmpe_nuisancechanged) replace nodraw
    matrix PMPE_CHANGED=r(plot_data)
    mata: st_numscalar("PMPE_error",all(st_matrix("PMPE_CHANGED"):==st_matrix("PMPE_NI")))
    test_assert "[current `family'] nuisance edits leave frozen influence values unchanged" scalar(PMPE_error)==1

}


capture mata: mata drop pmpe_influence_oracle()
_rdistnn_test_e_snapshot_drop, name(PMPE_ESTATE)
macro drop PMPE_*
