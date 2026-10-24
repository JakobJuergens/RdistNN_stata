// Adversarial postestimation audit. Run from the repository root after adding
// src/ to adopath. The parent audit serializes licensed Stata execution.
// Presence is checked through r(matrices), never assignment _rc: an absent
// r(matrix) can be assigned successfully as a missing 1-by-1 matrix in Stata.
version 17
do "tests/test_framework.do"
capture mkdir "build"
capture mkdir "build/premerge"

do "tests/helpers/postestimation-state.do"

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
    local panel_columns : colnames PMPE_PANELS
    local expectedpanels "estimates diagnostics weights_cumulative tuning"
    if "`family'"=="np" local expectedpanels "estimates diagnostics outcome weights_cumulative tuning"
    local expectedn : word count `expectedpanels'
    test_assert "[premerge `family'] dashboard proves plot-data absence and exact panel matrix" `hasdata'==0 & `haspanels'>0 & `nr'==. & `npanels'==`expectedn' & rowsof(PMPE_PANELS)==`expectedn' & colsof(PMPE_PANELS)==2 & "`panels'"=="`expectedpanels'"
    test_assert "[premerge `family'] dashboard counts are actual panel rows" PMPE_PANELS[1,2]==5 & PMPE_PANELS[`expectedn',2]==9 & "`panel_columns'"=="panel_id n_rows"

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

test_section "selected reporting weights only the requested targets"
do "tests/helpers/reporting-instrument.do"
capture program drop _pmpe_selected_weight_work
program define _pmpe_selected_weight_work
    tempfile root
    local shadow "`root'_sources"
    mkdir "`shadow'"
    capture noisily {
        // Count the columns actually supplied to the real weight primitives;
        // source shadows keep production counters and test hooks out of core.
        foreach module in model_mata plot_mata {
            mata: rdistnn_test_instrument("src/_rdistnn_"+st_local("module")+".ado",st_local("shadow")+"/step.ado","rdistnn_weights_from_rmat(R,","rdistnn_weights_from_rmat(rdistnn_test_weight_columns(R),")
            mata: rdistnn_test_instrument(st_local("shadow")+"/step.ado",st_local("shadow")+"/_rdistnn_"+st_local("module")+".ado","rdistnn_tdnn_weights(R,","rdistnn_tdnn_weights(rdistnn_test_weight_columns(R),")
        }
        mata: rdistnn_test_instrument("src/_rdistnn_model_decode.ado",st_local("shadow")+"/_rdistnn_model_decode.ado","X=st_matrix("+char(34)+"e(model_X)"+char(34)+"); T=st_matrix","st_numscalar("+char(34)+"rdistnn_test_model_validations"+char(34)+",st_numscalar("+char(34)+"rdistnn_test_model_validations"+char(34)+")+1); X=st_matrix("+char(34)+"e(model_X)"+char(34)+"); T=st_matrix")
        adopath ++ "`shadow'"
        scalar rdistnn_test_weighted_columns=0
        scalar rdistnn_test_model_validations=0
        clear
        set obs 12
        gen double pmpe_x=mod(_n,6)-3
        gen double pmpe_y=_n^2
        gen double pmpe_den=1+_n/20
        matrix PMPE_WEIGHT_TARGETS=(0\1\-.5)
        foreach family in np cate clate tdnn {
            local command "`family'_dnn"
            local opts "tuningmode(grid) sgrid(2 4 2) sdefault(4)"
            if "`family'"=="cate" local opts "`opts' pseudooutcome"
            if "`family'"=="clate" local opts "`opts' pseudooutcome denominator(pmpe_den)"
            if "`family'"=="tdnn" {
                local command "np_tdnn"
                local opts "s1(2) s2(4)"
            }
            quietly `command' pmpe_y, xcols(pmpe_x) x0matrix(PMPE_WEIGHT_TARGETS) `opts'
            foreach option in validated validated(1) {
                capture quietly `command', `option'
                local rejected=_rc
                test_assert "`family': replay exposes no private verification option" `rejected'==198
                capture quietly `command' pmpe_y, xcols(pmpe_x) x0matrix(PMPE_WEIGHT_TARGETS) `opts' `option'
                local rejected=_rc
                test_assert "`family': fresh estimation exposes no private verification option" `rejected'==198
                capture quietly rdistnn_predict, x0(.25) post `option'
                local rejected=_rc
                test_assert "`family': public posting exposes no private verification option" `rejected'==198
            }
            scalar rdistnn_test_model_validations=0
            quietly `command'
            test_assert "`family': ordinary replay validates its current snapshot" scalar(rdistnn_test_model_validations)==1
            quietly rdistnn_results, type(weights)
            matrix PMPE_FULL_WEIGHTS=r(data)
            _rdistnn_test_e_snapshot, name(PMPE_WEIGHT_STATE)
            scalar rdistnn_test_weighted_columns=0
            scalar rdistnn_test_model_validations=0
            quietly rdistnn_results, type(weights) targets(3 1 3)
            local weighted=scalar(rdistnn_test_weighted_columns)
            local validations=scalar(rdistnn_test_model_validations)
            matrix PMPE_SELECTED_WEIGHTS=r(data)
            test_assert "`family': selected export weights two columns" `weighted'==2
            test_assert "`family': selected export validates unchanged state once" `validations'==1
            mata: st_numscalar("PMPE_WEIGHT_EQUAL",all(st_matrix("PMPE_SELECTED_WEIGHTS"):==st_matrix("PMPE_FULL_WEIGHTS")[(25::36)\(1::12),]))
            test_assert "`family': selected weights retain original targets and exact values" scalar(PMPE_WEIGHT_EQUAL)
            _rdistnn_test_e_compare, snapshot(PMPE_WEIGHT_STATE)
            test_assert "`family': selected export preserves complete fitted state" r(equal)
            scalar rdistnn_test_model_validations=0
            quietly rdistnn_results, type(weights) targets(1)
            test_assert "`family': independent calls validate independently" scalar(rdistnn_test_model_validations)==1
            scalar rdistnn_test_weighted_columns=0
            scalar rdistnn_test_model_validations=0
            quietly `command'_plot, type(weights) targets(3 1) name(pmpe_selected) replace nodraw
            local weighted=scalar(rdistnn_test_weighted_columns)
            local validations=scalar(rdistnn_test_model_validations)
            test_assert "`family': real selected plot weights two columns" `weighted'==2
            test_assert "`family': real selected plot validates unchanged state once" `validations'==1
            _rdistnn_test_e_compare, snapshot(PMPE_WEIGHT_STATE)
            test_assert "`family': selected plot preserves complete fitted state" r(equal)
            scalar rdistnn_test_model_validations=0
            quietly `command'_plot, type(diagnostics) targets(3 1) name(pmpe_selected) replace nodraw
            test_assert "`family': diagnostic plots reuse current validation" scalar(rdistnn_test_model_validations)==1
            scalar rdistnn_test_model_validations=0
            quietly rdistnn_predict, x0(.25) sefit
            test_assert "`family': prediction validates source and newly posted state" scalar(rdistnn_test_model_validations)==2
            if "`family'"!="tdnn" {
                scalar rdistnn_test_model_validations=0
                quietly rdistnn_candidate, candidate(3)
                test_assert "`family': candidate validates source and newly posted state" scalar(rdistnn_test_model_validations)==2
            }
            scalar rdistnn_test_model_validations=0
            mata: PMPE_PRIVATE_WEIGHTS=rdistnn_model_weights((3,1))
            test_assert "`family': standalone model helper retains validation by default" scalar(rdistnn_test_model_validations)==1
            _rdistnn_test_e_snapshot_drop, name(PMPE_WEIGHT_STATE)
        }
    }
    local rc=_rc
    capture adopath - "`shadow'"
    quietly run "src/_rdistnn_model_mata.ado"
    quietly run "src/_rdistnn_plot_mata.ado"
    foreach file in step _rdistnn_model_mata _rdistnn_plot_mata _rdistnn_model_decode {
        capture erase "`shadow'/`file'.ado"
    }
    rmdir "`shadow'"
    capture scalar drop rdistnn_test_weighted_columns
    capture scalar drop rdistnn_test_model_validations
    _rdistnn_test_e_snapshot_drop, name(PMPE_WEIGHT_STATE)
    if `rc' exit `rc'
end
_pmpe_selected_weight_work
program drop _pmpe_selected_weight_work

test_section "all-target reporting follows matrix capacity, not numlist capacity"
mata: st_numscalar("PMPE_NUMERIC_ALL",all(rdistnn_model_target_index(J(1,0,.),2501):==(1..2501)))
test_assert "numeric all-target selection exceeds the numlist limit" scalar(PMPE_NUMERIC_ALL)
mata: st_numscalar("PMPE_NUMERIC_SELECTED",all(rdistnn_model_target_index((3,1),5):==(3,1)))
test_assert "numeric target selection preserves requested order" scalar(PMPE_NUMERIC_SELECTED)
// BE cannot represent this fitted result. SE/MP exercise the real public path
// once; small selected-target/plot contracts above cover every shared consumer.
if c(max_matsize)>=2501 {
    local pmpe_matsize=c(matsize)
    if c(matsize)<2501 set matsize 2501
    clear
    set obs 8
    gen double pmpe_large_x=_n
    gen double pmpe_large_y=_n^2
    mata: st_matrix("PMPE_LARGE_TARGETS",100:+(1::2501):/1000)
    quietly np_dnn pmpe_large_y, xcols(pmpe_large_x) ///
        x0matrix(PMPE_LARGE_TARGETS) s(3) nose
    capture quietly rdistnn_results
    local pmpe_large_rc=_rc
    test_assert "valid large fit reports all targets" `pmpe_large_rc'==0
    if !`pmpe_large_rc' {
        matrix PMPE_LARGE_INDEX=r(index)
        matrix PMPE_LARGE_FIT=r(fit)
        mata: st_numscalar("PMPE_LARGE_EQUAL",all(st_matrix("PMPE_LARGE_FIT"):==st_matrix("e(estimates)")) & all(st_matrix("PMPE_LARGE_INDEX")[,1]:==(1::2501)))
        test_assert "large report preserves all target identities and estimates" scalar(PMPE_LARGE_EQUAL) & rowsof(PMPE_LARGE_INDEX)==2501
    }
    quietly rdistnn_results, targets(2501 1 2501)
    matrix PMPE_LARGE_SELECTED=r(index)
    test_assert "explicit large-target IDs retain ordered deduplication" rowsof(PMPE_LARGE_SELECTED)==2 & PMPE_LARGE_SELECTED[1,1]==2501 & PMPE_LARGE_SELECTED[2,1]==1
    test_assert "large reporting preserves the fit and caller data" e(m)==2501 & "`e(cmd)'"=="np_dnn" & _N==8 & pmpe_large_y[8]==64
    set matsize `pmpe_matsize'
}

// Remove audit-only programs/global snapshots; production e() remains active.
capture program drop _pmpe_snapshot
capture program drop _pmpe_unchanged
_rdistnn_test_e_snapshot_drop, name(PMPE_ESTATE)
macro drop PMPE_*
