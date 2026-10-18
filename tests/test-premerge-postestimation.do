// Adversarial postestimation audit. Run from the repository root after adding
// src/ to adopath. The parent audit serializes licensed Stata execution.
// Presence is checked through r(matrices), never assignment _rc: an absent
// r(matrix) can be assigned successfully as a missing 1-by-1 matrix in Stata.
version 17
do "tests/test_framework.do"
capture mkdir "build"
capture mkdir "build/premerge"

capture program drop _pmpe_snapshot
program define _pmpe_snapshot
    // Snapshot every stored matrix (including stripes), scalar and macro.
    // This creates no estimation stores or dataset sample-marker variables.
    global PMPE_matrices : e(matrices)
    global PMPE_scalars : e(scalars)
    global PMPE_macros : e(macros)
    global PMPE_datalabel : data label
    unab variables : _all
    global PMPE_variables "`variables'"
    local i 0
    foreach variable of local variables {
        local ++i
        global PMPE_type`i' : type `variable'
        global PMPE_format`i' : format `variable'
        global PMPE_label`i' : variable label `variable'
        global PMPE_valuelabel`i' : value label `variable'
    }
    local i 0
    foreach key of global PMPE_matrices {
        local ++i
        matrix PMPE_m`i' = e(`key')
    }
    local i 0
    foreach key of global PMPE_scalars {
        local ++i
        scalar PMPE_s`i' = e(`key')
    }
    local i 0
    foreach key of global PMPE_macros {
        local ++i
        global PMPE_g`i' `"`e(`key')'"'
    }
    tempvar sample
    quietly gen byte `sample' = e(sample)
    quietly mkmat `sample', matrix(PMPE_ESAMPLE)
end

capture program drop _pmpe_unchanged
program define _pmpe_unchanged, rclass
    syntax, DATA(string)
    local ok 1
    capture cf _all using `"`data'"'
    if _rc local ok 0
    local datalabel : data label
    if `"`datalabel'"'!=`"$PMPE_datalabel"' local ok 0
    unab variables : _all
    if "`variables'"!="$PMPE_variables" local ok 0
    local i 0
    foreach variable of local variables {
        local ++i
        foreach property in type format label valuelabel {
            local original : copy global PMPE_`property'`i'
            if "`property'"=="type" local current : type `variable'
            if "`property'"=="format" local current : format `variable'
            if "`property'"=="label" local current : variable label `variable'
            if "`property'"=="valuelabel" local current : value label `variable'
            if `"`current'"'!=`"`original'"' local ok 0
        }
    }
    local matrices : e(matrices)
    local scalars : e(scalars)
    local macros : e(macros)
    if "`matrices'"!="$PMPE_matrices" local ok 0
    if "`scalars'"!="$PMPE_scalars" local ok 0
    if "`macros'"!="$PMPE_macros" local ok 0
    local i 0
    foreach key of global PMPE_matrices {
        local ++i
        tempname current
        matrix `current' = e(`key')
        if rowsof(`current')!=rowsof(PMPE_m`i') | colsof(`current')!=colsof(PMPE_m`i') local ok 0
        else {
            mata: st_numscalar("PMPE_match", all(st_matrix("`current'"):==st_matrix("PMPE_m`i'")) & all(st_matrixrowstripe("`current'"):==st_matrixrowstripe("PMPE_m`i'")) & all(st_matrixcolstripe("`current'"):==st_matrixcolstripe("PMPE_m`i'")))
            if scalar(PMPE_match)!=1 local ok 0
        }
    }
    local i 0
    foreach key of global PMPE_scalars {
        local ++i
        if e(`key')!=scalar(PMPE_s`i') local ok 0
    }
    local i 0
    foreach key of global PMPE_macros {
        local ++i
        local original : copy global PMPE_g`i'
        if `"`e(`key')'"'!=`"`original'"' local ok 0
    }
    tempvar sample
    tempname current_sample
    quietly gen byte `sample' = e(sample)
    quietly mkmat `sample', matrix(`current_sample')
    if rowsof(`current_sample')!=rowsof(PMPE_ESAMPLE) local ok 0
    else {
        mata: st_numscalar("PMPE_match",all(st_matrix("`current_sample'"):==st_matrix("PMPE_ESAMPLE")))
        if scalar(PMPE_match)!=1 local ok 0
    }
    return scalar ok = `ok'
end

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

test_section "premerge postestimation: selectors, genuine returns and state"
foreach family in cate np clate {
    rdnn_dgp_basic, n(48) seed(923)
    gen double denominator = .8+.1*sin(x1)
    // Deliberate collisions with plot fields, plus a string and data label:
    // preserve/restore must leave these caller-owned columns in place.
    gen long target_id = _n
    gen double estimate = 9+_n
    gen str12 contribution = "caller"
    label data "premerge postestimation caller data"
    matrix PMPE_X0 = (.7,0 \ -.8,0 \ .1,0 \ 1.2,0 \ -.2,0)
    local opts ""
    if "`family'"=="cate" local opts "pseudooutcome"
    if "`family'"=="clate" local opts "pseudooutcome denominator(denominator)"
    quietly `family'_dnn y, xcols(x1 x2) x0matrix(PMPE_X0) ///
        tuningmode(grid) sgrid(3 7 3 48) sdefault(7) `opts'
    matrix PMPE_E = e(estimates)
    matrix PMPE_SE = e(se_targets)
    matrix PMPE_V = e(estimate_vcov)
    matrix PMPE_META = e(tuning_candidates)
    matrix PMPE_GE = e(tuning_estimates)
    matrix PMPE_GSE = e(tuning_se)
    matrix PMPE_DIAG = e(diagnostics)
    tempfile caller diagnostics_saved
    quietly save "`caller'", replace
    _pmpe_snapshot

    quietly `family'_dnn_diagnostics, targets(5 2 5) ///
        saving("`diagnostics_saved'") replace
    local matrices : r(matrices)
    local required "localization diagnostics compact"
    local present : list required in matrices
    local nr = r(n_rows)
    local selected "`r(selected)'"
    matrix PMPE_D = r(diagnostics)
    matrix PMPE_L = r(localization)
    if "`family'"=="clate" {
        local required "components status"
        local c_present : list required in matrices
        matrix PMPE_C = r(components)
        matrix PMPE_S = r(status)
        test_assert "[premerge `family'] component/status matrices actually returned with duplicate rows" `c_present' & rowsof(PMPE_C)==3 & colsof(PMPE_C)==9 & rowsof(PMPE_S)==3 & colsof(PMPE_S)==5 & PMPE_C[1,1]==5 & PMPE_C[2,1]==2 & PMPE_C[3,1]==5
    }
    local cols = cond("`family'"=="np",20,28)
    test_assert "[premerge `family'] diagnostics preserve duplicate IDs and requested order" `present' & `nr'==3 & "`selected'"=="5 2 5" & rowsof(PMPE_D)==3 & colsof(PMPE_D)==`cols' & rowsof(PMPE_L)==3 & colsof(PMPE_L)==14 & PMPE_D[1,1]==5 & PMPE_D[2,1]==2 & PMPE_D[3,1]==5
    mata: st_numscalar("PMPE_error", max(abs(st_matrix("PMPE_D")[,4..cols(st_matrix("PMPE_D"))]-st_matrix("PMPE_DIAG")[(5,2,5),])))
    test_assert "[premerge `family'] full diagnostic values follow original rows, including duplicates" scalar(PMPE_error)<1e-12
    preserve
    use "`diagnostics_saved'", clear
    local savecols = cond("`family'"=="clate",41,`cols')
    quietly describe
    local saveok = _N==3 & r(k)==`savecols' & target_id[1]==5 & target_id[2]==2 & target_id[3]==5
    unab variables : _all
    foreach variable of local variables {
        local storage : type `variable'
        if "`storage'"!="double" local saveok 0
    }
    mkmat _all, matrix(PMPE_SAVED)
    if "`family'"=="clate" matrix PMPE_EXPECTED = PMPE_D, PMPE_C[1...,2...], PMPE_S
    else matrix PMPE_EXPECTED = PMPE_D
    mata: st_numscalar("PMPE_error", all(st_matrix("PMPE_SAVED"):==st_matrix("PMPE_EXPECTED")))
    local saveok = `saveok' & scalar(PMPE_error)==1
    restore
    test_assert "[premerge `family'] saved file exactly matches selected numeric return, with duplicate identities" `saveok'

    quietly `family'_dnn_diagnostics, candidates(3 1 3)
    local matrices : r(matrices)
    local absent "localization components status"
    local leaked : list matrices & absent
    local selected "`r(selected)'"
    matrix PMPE_CD = r(diagnostics)
    test_assert "[premerge `family'] candidate scope has no target-only matrices and preserves duplicates" "`leaked'"=="" & "`selected'"=="3 1 3" & rowsof(PMPE_CD)==3 & colsof(PMPE_CD)==13 & PMPE_CD[1,1]==3 & PMPE_CD[2,1]==1 & PMPE_CD[3,1]==3
    quietly `family'_dnn_diagnostics, allcandidates
    matrix PMPE_ALLD = r(diagnostics)
    test_assert "[premerge `family'] infeasible and aliased candidate rows remain in diagnostics" rowsof(PMPE_ALLD)==4 & PMPE_ALLD[4,6]==0 & PMPE_ALLD[1,2]==PMPE_ALLD[3,2] & PMPE_ALLD[4,8]==48 & missing(PMPE_ALLD[4,9])

    quietly `family'_dnn, level(80)
    local matrices : r(matrices)
    local required "table summary"
    local present : list required in matrices
    matrix PMPE_T = r(table)
    matrix PMPE_SUM = r(summary)
    test_assert "[premerge `family'] replay actually returns full report matrices" `present' & rowsof(PMPE_T)==9 & colsof(PMPE_T)==5 & rowsof(PMPE_SUM)==5 & colsof(PMPE_SUM)==6
    quietly `family'_dnn_plot, targets(5 2 5) level(80) name(pmpe_est) replace nodraw
    local matrices : r(matrices)
    local hasdata : list posof "plot_data" in matrices
    local nr = r(n_rows)
    local targets "`r(targets)'"
    local selected "`r(selected_targets)'"
    matrix PMPE_P = r(plot_data)
    test_assert "[premerge `family'] estimate plot deduplicates and filters in stored target order" `hasdata'>0 & `nr'==2 & rowsof(PMPE_P)==2 & colsof(PMPE_P)==10 & PMPE_P[1,6]==2 & PMPE_P[2,6]==5 & "`targets'"=="2 5" & "`selected'"=="2 5"
    test_assert "[premerge `family'] replay and plot bounds match exactly at requested level" abs(PMPE_P[1,9]-PMPE_T[5,2])<1e-12 & abs(PMPE_P[2,10]-PMPE_T[6,5])<1e-12 & abs(PMPE_P[1,8]^2-PMPE_V[2,2])<1e-12
    quietly `family'_dnn_plot, targets(5 2 5) candidates(3 1 3) ///
        level(80) name(pmpe_grid) replace nodraw
    local nr = r(n_rows)
    local targets "`r(targets)'"
    local candidates "`r(candidates)'"
    matrix PMPE_GP = r(plot_data)
    test_assert "[premerge `family'] candidate estimates deduplicate with requested candidate/stored target order" `nr'==4 & rowsof(PMPE_GP)==4 & colsof(PMPE_GP)==10 & PMPE_GP[1,3]==3 & PMPE_GP[3,3]==1 & PMPE_GP[1,6]==2 & PMPE_GP[2,6]==5 & "`targets'"=="2 5" & "`candidates'"=="1 3"
    test_assert "[premerge `family'] aliased candidate estimates and SEs agree" PMPE_GP[1,7]==PMPE_GP[3,7] & PMPE_GP[2,8]==PMPE_GP[4,8] & abs(PMPE_GP[1,7]-PMPE_GE[3,2])<1e-12 & abs(PMPE_GP[1,9]-(PMPE_GE[3,2]-invnormal(.9)*PMPE_GSE[3,2]))<1e-12
    quietly `family'_dnn_plot, type(tuning) targets(5 2 5) candidates(3 1 3) ///
        name(pmpe_tune) replace nodraw
    local nr = r(n_rows)
    local selected "`r(selected_targets)'"
    matrix PMPE_TP = r(plot_data)
    test_assert "[premerge `family'] tuning loops preserve requested unique target/candidate order" `nr'==4 & colsof(PMPE_TP)==9 & PMPE_TP[1,1]==3 & PMPE_TP[3,1]==1 & PMPE_TP[1,5]==5 & PMPE_TP[2,5]==2 & "`selected'"=="5 2"
    quietly `family'_dnn_plot, type(diagnostics) targets(5 2 5) name(pmpe_diag) replace nodraw
    local nr = r(n_rows)
    matrix PMPE_DP = r(plot_data)
    local dcols = cond("`family'"=="np",7,9)
    test_assert "[premerge `family'] diagnostic plot preserves requested unique order and canonical identity" `nr'==2 & rowsof(PMPE_DP)==2 & colsof(PMPE_DP)==`dcols' & PMPE_DP[1,1]==5 & PMPE_DP[2,1]==2 & PMPE_DP[1,`=`dcols'-1']==5 & PMPE_DP[1,`dcols']==0

    quietly `family'_dnn_plot, style(profile) xindex(1) level(80) ///
        name(pmpe_export) replace nodraw ///
        saving("build/premerge/postestimation-`family'-interval-profile.svg")
    local graphtype "`r(type)'"
    local graphview "`r(view)'"
    local graphname "`r(graph)'"
    local graphsaved `"`r(saved)'"'
    matrix PMPE_EXPORT = r(plot_data)
    capture confirm file "build/premerge/postestimation-`family'-interval-profile.svg"
    local export_exists = _rc==0
    test_assert "[premerge `family'] profile export carries graph/type/view/saved metadata" `export_exists' & "`graphtype'"=="estimates" & "`graphview'"=="profile" & "`graphname'"=="pmpe_export" & `"`graphsaved'"'=="build/premerge/postestimation-`family'-interval-profile.svg"
    preserve
    clear
    svmat double PMPE_EXPORT, names(col)
    export delimited using "build/premerge/postestimation-`family'-interval-profile-data.csv", replace
    restore
    if "`family'"=="np" {
        quietly np_dnn_plot, type(outcome) xindex(1) name(pmpe_outcome_export) replace nodraw ///
            saving("build/premerge/postestimation-np-outcome-profile.svg")
        matrix PMPE_OUTCOME = r(plot_data)
        test_assert "[premerge NP] outcome export retains exact mean and quantile geometry" rowsof(PMPE_OUTCOME)==5 & colsof(PMPE_OUTCOME)==8 & abs(PMPE_OUTCOME[3,4]-PMPE_E[1,3])<1e-12 & PMPE_OUTCOME[3,6]<=PMPE_OUTCOME[3,7] & PMPE_OUTCOME[3,7]<=PMPE_OUTCOME[3,8]
        preserve
        clear
        svmat double PMPE_OUTCOME, names(col)
        export delimited using "build/premerge/postestimation-np-outcome-profile-data.csv", replace
        restore
    }

    local original_level = c(level)
    set level 83
    quietly `family'_dnn
    local default_level = r(level)
    matrix PMPE_DEFAULT_REPORT = r(table)
    quietly `family'_dnn_plot, name(pmpe_default_level) replace nodraw
    matrix PMPE_DEFAULT_PLOT = r(plot_data)
    set level `original_level'
    test_assert "[premerge `family'] configured c(level) is shared by replay and plots" `default_level'==83 & abs(PMPE_DEFAULT_REPORT[8,1]-invnormal(.915))<1e-12 & abs(PMPE_DEFAULT_PLOT[1,9]-PMPE_DEFAULT_REPORT[5,1])<1e-12 & abs(PMPE_DEFAULT_PLOT[1,10]-PMPE_DEFAULT_REPORT[6,1])<1e-12

    quietly `family'_dnn_plot, type(dashboard) name(pmpe_dash) replace nodraw
    local matrices : r(matrices)
    local hasdata : list posof "plot_data" in matrices
    local haspanels : list posof "panel_info" in matrices
    local nr = r(n_rows)
    local npanels = r(n_panels)
    local panels "`r(panels)'"
    local panels : list retokenize panels
    matrix PMPE_PANELS = r(panel_info)
    local expectedpanels "estimates diagnostics weights_cumulative tuning"
    if "`family'"=="np" local expectedpanels "estimates diagnostics outcome weights_cumulative tuning"
    local expectedn : word count `expectedpanels'
    test_assert "[premerge `family'] dashboard proves plot-data absence and exact panel matrix" `hasdata'==0 & `haspanels'>0 & `nr'==. & `npanels'==`expectedn' & rowsof(PMPE_PANELS)==`expectedn' & colsof(PMPE_PANELS)==3 & "`panels'"=="`expectedpanels'"
    test_assert "[premerge `family'] dashboard counts reflect rows vs legacy counts" PMPE_PANELS[1,2]==5 & PMPE_PANELS[1,3]==5 & PMPE_PANELS[`expectedn',2]==9 & PMPE_PANELS[`expectedn',3]==3

    capture noisily `family'_dnn_diagnostics, targets(2) saving("build/premerge/postestimation-no-such-dir/report.dta")
    local save_rc = _rc
    test_assert "[premerge `family'] failed diagnostic save reports failure" `save_rc'!=0
    capture noisily `family'_dnn_plot, saving("build/premerge/postestimation-no-such-dir/plot.svg") ///
        name(pmpe_badexport) replace nodraw
    local export_rc = _rc
    test_assert "[premerge `family'] failed graph export reports failure" `export_rc'!=0
    capture noisily `family'_dnn_plot, type(weights) candidates(4) name(pmpe_badcandidate) replace nodraw
    local candidate_rc = _rc
    test_assert "[premerge `family'] explicit infeasible plot candidate rejected" `candidate_rc'==498
    _pmpe_unchanged, data("`caller'")
    local unchanged = r(ok)
    test_assert "[premerge `family'] all successes and I/O failures preserve all data and every e() value/stripe" `unchanged'==1
}

test_section "premerge postestimation: original sample geometry and score reconstruction"
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
    test_assert "[premerge `family'] kernel distances/ranks/weights/cumulative shares match independent original-row oracle" scalar(PMPE_geometry_ok)==1
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
    // Same count and unchanged X are insufficient if an input signal changes.
    replace y = y+.01 in 20
    capture noisily `family'_dnn_plot, type(weights) view(kernel) name(pmpe_changed) replace nodraw
    local change_rc = _rc
    test_assert "[premerge `family'] in-sample signal change rejected even for geometric observation view" `change_rc'==498
    quietly use "`restricted'", clear
    // Re-loading the original data keeps the current e(sample) uncertain on
    // some Stata versions; same-session sorting is tested using a fresh fit.
    quietly `family'_dnn y if source_id>=5 & source_id<=43 in 3/46, ///
        xcols(x1 x2) x0matrix(PMPE_X0) s(6) `opts'
    gsort -source_id
    capture noisily `family'_dnn_plot, type(weights) view(kernel) name(pmpe_order) replace nodraw
    local order_rc = _rc
    test_assert "[premerge `family'] sorting that carries the sample marker still fails the order signature" `order_rc'==498
    quietly `family'_dnn_plot, type(estimates) name(pmpe_stored) replace nodraw
    local matrices : r(matrices)
    local present : list posof "plot_data" in matrices
    matrix PMPE_STORED = r(plot_data)
    test_assert "[premerge `family'] stored estimate view remains independent of changed dataset order" `present'>0 & rowsof(PMPE_STORED)==2 & colsof(PMPE_STORED)==10
}

test_section "premerge postestimation: supplied nuisances are guarded and reconstructed"
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
    scalar PMPE_original_m0 = m0[1]
    replace m0 = m0+.001 in 1
    capture noisily `family'_dnn_plot, type(influence) name(pmpe_nuisancebad) replace nodraw
    local nuisance_rc = _rc
    test_assert "[premerge `family'] changing original nuisance is rejected" `nuisance_rc'==498
    replace m0 = scalar(PMPE_original_m0) in 1
    capture noisily `family'_dnn_plot, type(influence) name(pmpe_nuisanceback) replace nodraw
    local restored_rc = _rc
    test_assert "[premerge `family'] restoring exact nuisance bytes re-enables observation views" `restored_rc'==0
}

test_section "premerge postestimation: cold restored summaries and independent plot axes"
foreach family in cate np clate {
    rdnn_dgp_basic, n(32) seed(901)
    rename x1 x0_1
    rename x2 target_id
    gen double denominator = 1+.05*sin(x0_1)
    matrix PMPE_COLLISION = (1,0 \ -.7,0 \ .4,0 \ -.2,0)
    local opts ""
    if "`family'"=="cate" local opts "pseudooutcome"
    if "`family'"=="clate" local opts "pseudooutcome denominator(denominator)"
    quietly `family'_dnn y, xcols(x0_1 target_id) x0matrix(PMPE_COLLISION) ///
        s(5) `opts' nose
    tempfile cold
    quietly estimates save "`cold'", replace
    matrix PMPE_BEFORE = e(estimates)
    discard
    quietly estimates use "`cold'"
    quietly `family'_dnn, level(90)
    local matrices : r(matrices)
    local required "table summary"
    local present : list required in matrices
    matrix PMPE_COLD = r(table)
    test_assert "[premerge `family'] cold saved current fit returns full nose report" `present' & rowsof(PMPE_COLD)==9 & colsof(PMPE_COLD)==4 & PMPE_COLD[1,3]==PMPE_BEFORE[1,3] & missing(PMPE_COLD[2,3]) & missing(PMPE_COLD[5,3]) & missing(PMPE_COLD[6,3])
    quietly `family'_dnn_plot, xindex(1) style(profile) name(pmpe_cold) replace nodraw ///
        saving("build/premerge/postestimation-`family'-nose-profile.svg")
    local matrices : r(matrices)
    local present : list posof "plot_data" in matrices
    local view "`r(view)'"
    local saved `"`r(saved)'"'
    matrix PMPE_CP = r(plot_data)
    test_assert "[premerge `family'] collision names preserve coordinates and canonical target IDs" `present'>0 & rowsof(PMPE_CP)==4 & colsof(PMPE_CP)==10 & PMPE_CP[1,1]==1 & PMPE_CP[2,1]==-.7 & PMPE_CP[1,6]==1 & PMPE_CP[2,6]==2 & "`view'"=="profile" & missing(PMPE_CP[1,8]) & missing(PMPE_CP[1,9]) & missing(PMPE_CP[1,10])
    capture confirm file "build/premerge/postestimation-`family'-nose-profile.svg"
    local exists = _rc==0
    test_assert "[premerge `family'] cold-session SVG export has exact saved metadata" `exists' & `"`saved'"'=="build/premerge/postestimation-`family'-nose-profile.svg"
    quietly `family'_dnn_diagnostics, targets(4 2)
    local matrices : r(matrices)
    local present : list posof "localization" in matrices
    matrix PMPE_CL = r(localization)
    test_assert "[premerge `family'] cold restored diagnostics retain supplied coordinate semantics" `present'>0 & rowsof(PMPE_CL)==2 & colsof(PMPE_CL)==14 & PMPE_CL[1,1]==4 & PMPE_CL[1,2]==-.2 & PMPE_CL[2,1]==2 & PMPE_CL[2,2]==-.7
}

test_section "premerge postestimation: target-specific ratio inference masks"
clear
set obs 6
gen double x = _n
gen double numerator = cond(_n==6,4,1)
gen double denominator = 1
matrix PMPE_MIXED = (-1 \ 8)
quietly clate_dnn numerator, denominator(denominator) pseudooutcome ///
    xcols(x) x0matrix(PMPE_MIXED) tuningmode(grid) sgrid(5 6) sdefault(5)
matrix PMPE_MIXED_REPORT = r(table)
matrix PMPE_MIXED_STATUS = e(status)
matrix PMPE_MIXED_V = e(estimate_vcov)
test_assert "[premerge CLATE] mixed targets retain zero covariance separately from positive Wald variance" PMPE_MIXED_STATUS[1,4]==2 & PMPE_MIXED_STATUS[1,5]==0 & PMPE_MIXED_STATUS[2,5]==1 & PMPE_MIXED_V[1,1]==0 & PMPE_MIXED_V[1,2]==0 & PMPE_MIXED_V[2,2]>0 & missing(PMPE_MIXED_REPORT[5,1]) & !missing(PMPE_MIXED_REPORT[5,2])
foreach kind in estimates tuning {
    local candidate_option ""
    if "`kind'"=="tuning" local candidate_option "candidates(1)"
    quietly clate_dnn_plot, type(`kind') level(95) `candidate_option' ///
        name(pmpe_mixed) replace nodraw
    matrix PMPE_MIXED_PLOT = r(plot_data)
    local columns : colnames PMPE_MIXED_PLOT
    local sec : list posof "se" in columns
    local lowerc : list posof "ci_lower" in columns
    local upperc : list posof "ci_upper" in columns
    test_assert "[premerge CLATE] `kind' masks only the zero-variance target, retaining its zero SE" rowsof(PMPE_MIXED_PLOT)==2 & PMPE_MIXED_PLOT[1,`sec']==0 & missing(PMPE_MIXED_PLOT[1,`lowerc']) & missing(PMPE_MIXED_PLOT[1,`upperc']) & abs(PMPE_MIXED_PLOT[2,`lowerc']-PMPE_MIXED_REPORT[5,2])<1e-12 & abs(PMPE_MIXED_PLOT[2,`upperc']-PMPE_MIXED_REPORT[6,2])<1e-12
}

test_section "premerge postestimation: undefined selection has no phantom matrix"
clear
set obs 20
gen double x = _n
gen double numerator = 1+sin(_n)
gen double denominator = 0
matrix PMPE_UNDEFINED = (0 \ 10 \ 21)
quietly clate_dnn numerator, denominator(denominator) pseudooutcome ///
    xcols(x) x0matrix(PMPE_UNDEFINED) tuningmode(grid)
foreach kind in estimates tuning {
    quietly clate_dnn_plot, type(`kind') targets(3 1 3) ///
        name(pmpe_undefined) replace nodraw
    local matrices : r(matrices)
    local hasdata : list posof "plot_data" in matrices
    local nr = r(n_rows)
    local omitted = r(n_omitted)
    local selected "`r(selected_targets)'"
    local expected = cond("`kind'"=="estimates",2,2*e(tuning_n_feasible))
    local expectedselected "1 3"
    if "`kind'"=="tuning" local expectedselected "3 1"
    test_assert "[premerge CLATE] wholly undefined `kind' selection returns zero rows and no plot_data member" `hasdata'==0 & `nr'==0 & `omitted'==`expected' & "`selected'"=="`expectedselected'"
}

// Remove audit-only programs/global snapshots; production e() remains active.
capture program drop _pmpe_snapshot
capture program drop _pmpe_unchanged
capture mata: mata drop pmpe_influence_oracle()
macro drop PMPE_*
