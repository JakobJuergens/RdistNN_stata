// Common plot contract: selected data, scientific bounds, validation and guards.
version 17
do "tests/test_framework.do"
quietly do "tests/fixtures/r_reporting_fixtures.do"
capture mkdir "build"
capture mkdir "build/coherence"
test_section "estimator plot coherence: selected data and confidence bounds"
foreach family in cate np clate {
    rdnn_dgp_basic, n(60) seed(521)
    gen double sy = y
    gen double sd = .8 + .1*x1
    matrix PX = (-1,0 \ -.5,0 \ 0,0 \ .5,0 \ 1,0)
    if "`family'"=="cate" quietly cate_dnn y, xcols(x1 x2) x0matrix(PX) tcol(D) ecol(e) m0col(m0) m1col(m1) tuningmode(grid)
    if "`family'"=="np" quietly np_dnn y, xcols(x1 x2) x0matrix(PX) tuningmode(grid)
    if "`family'"=="clate" quietly clate_dnn sy, xcols(x1 x2) x0matrix(PX) denominator(sd) pseudooutcome tuningmode(grid)
    estimates store plot_`family'
    matrix PB = e(b)
    local original_n = _N
    capture noisily `family'_dnn_plot, targets(2 4) level(90) name(plot_contract) replace nodraw
    local rc = _rc
    local nr = r(n_rows)
    capture matrix PD = r(plot_data)
    local drc = _rc
    test_assert "[plots `family'] selected estimate data exposed" `rc'==0 & `drc'==0 & `nr'==2
    if `drc'==0 {
        test_assert "[plots `family'] selected target identities retained" PD[1,6]==2 & PD[2,6]==4
        test_assert "[plots `family'] requested confidence level used" abs(PD[1,9]-(PD[1,7]-invnormal(.95)*PD[1,8]))<1e-12 & PD[1,10]>PD[1,9]
    }
    capture noisily `family'_dnn_plot, type(tuning) targets(2 4) candidates(1 3) name(plot_contract) replace nodraw
    local nr = r(n_rows)
    local legacy_n = r(n_plotted)
    capture matrix PT = r(plot_data)
    local drc = _rc
    test_assert "[plots `family'] tuning rows and legacy candidate count distinguished" `drc'==0 & `nr'==4 & `legacy_n'==2
    if `drc'==0 test_assert "[plots `family'] tuning identities retained" PT[1,1]==1 & PT[3,1]==3 & PT[1,5]==2 & PT[2,5]==4
    capture noisily `family'_dnn_plot, type(diagnostics) targets(2 4) name(plot_contract) replace nodraw
    local nr = r(n_rows)
    capture matrix PD = r(plot_data)
    local drc = _rc
    test_assert "[plots `family'] selected diagnostics returned" `drc'==0 & `nr'==2
    if `drc'==0 test_assert "[plots `family'] diagnostics target identity" PD[1,1]==2 & PD[2,1]==4
    capture noisily `family'_dnn_plot, type(diagnostics) candidates(1 3) name(plot_contract) replace nodraw
    local nr = r(n_rows)
    capture matrix PD = r(plot_data)
    local drc = _rc
    test_assert "[plots `family'] candidate diagnostics returned" `drc'==0 & `nr'==2
    if `drc'==0 {
        local cn : colnames PD
        local ci : list posof "candidate_id" in cn
        test_assert "[plots `family'] diagnostics canonical candidate identity" `ci'>0
        if `ci'>0 test_assert "[plots `family'] candidate rows preserve grid ids" PD[1,`ci']==1 & PD[2,`ci']==3
    }
    capture noisily `family'_dnn_plot, type(weights) view(profile) candidates(1 3) logweights name(plot_contract) replace nodraw
    local nr = r(n_rows)
    capture matrix PW = r(plot_data)
    local drc = _rc
    test_assert "[plots `family'] rank weights expose plotted positive rows" `drc'==0 & `nr'>0 & `nr'<600
    if `drc'==0 & colsof(PW)>=2 {
        mata: st_numscalar("plot_min_weight", min(st_matrix("PW")[,2]))
        test_assert "[plots `family'] log rank data exclude nonpositive rows" scalar(plot_min_weight)>0 & rowsof(PW)==`nr'
    }
    capture noisily `family'_dnn_plot, type(weights) view(kernel) targets(2 4) name(plot_contract) replace nodraw
    local nr = r(n_rows)
    capture matrix PK = r(plot_data)
    local drc = _rc
    test_assert "[plots `family'] kernel rows exposed" `drc'==0 & `nr'==120
    if `drc'==0 test_assert "[plots `family'] kernel target identity" PK[1,1]==2 & PK[61,1]==4
    capture noisily `family'_dnn_plot, type(weights) targets(2) replace nodraw
    test_assert "[plots `family'] rank target selection is supported" _rc==0
    matrix PRSEL=r(plot_data)
    test_assert "[plots `family'] rank profile exposes target-specific rows" rowsof(PRSEL)==60 & colsof(PRSEL)==7
    capture noisily `family'_dnn_plot, type(diagnostics) targets(2) candidates(1) replace nodraw
    test_assert "[plots `family'] rejects conflicting diagnostics selection" _rc==198
    capture noisily `family'_dnn_plot, type(tuning) style(forest) replace nodraw
    test_assert "[plots `family'] rejects unused style option" _rc==198
    capture noisily `family'_dnn_plot, type(weights) view(cumulative) logweights replace nodraw
    test_assert "[plots `family'] rejects unused logweights option" _rc==198
    capture noisily `family'_dnn_plot, type(dashboard) targets(2) replace nodraw
    test_assert "[plots `family'] rejects unused dashboard selection" _rc==198
    capture noisily `family'_dnn_plot, type(dashboard) name(plot_contract) replace nodraw
    local nr = r(n_rows)
    local npanels = r(n_panels)
    local panels "`r(panels)'"
    local returned : r(matrices)
    local has_plot_data : list posof "plot_data" in returned
    local has_panel_info : list posof "panel_info" in returned
    capture matrix PAN = r(panel_info)
    local prc = _rc
    test_assert "[plots `family'] dashboard returns panel metadata without single data table" `has_panel_info'>0 & `has_plot_data'==0 & `prc'==0 & colsof(PAN)>=3 & `npanels'>=4 & `nr'==. & "`panels'"!=""
    quietly `family'_dnn_plot, style(profile) name(plot_contract) replace nodraw saving("build/coherence/`family'-profile.svg")
    quietly `family'_dnn_plot, style(profile) name(plot_contract) replace nodraw saving("build/coherence/`family'-profile.png")
    capture confirm file "build/coherence/`family'-profile.svg"
    test_assert "[plots `family'] actual SVG export" _rc==0
    if "`family'"=="np" {
        quietly np_dnn_plot, type(outcome) targets(2 4) name(plot_contract) replace nodraw
        local nr = r(n_rows)
        capture matrix PO = r(plot_data)
        local drc = _rc
        test_assert "[plots np] selected outcome rows exposed" `drc'==0 & `nr'==2
        if `drc'==0 test_assert "[plots np] outcome target ids retained" PO[1,3]==2 & PO[2,3]==4
        quietly np_dnn_plot, type(outcome) name(plot_contract) replace nodraw saving("build/coherence/np-outcome.svg")
        quietly np_dnn_plot, type(outcome) name(plot_contract) replace nodraw saving("build/coherence/np-outcome.png")
    }
    else {
        if "`family'"=="cate" {
        quietly `family'_dnn_plot, type(overlap) targets(2 4) name(plot_contract) replace nodraw
        local nr = r(n_rows)
        capture matrix PO = r(plot_data)
        local drc = _rc
        // Paired-signal CLATE has no instrument overlap; use the supplied-nuisance route below.
        test_assert "[plots `family'] selected overlap rows exposed" `drc'==0 & `nr'==2
        }
        quietly `family'_dnn_plot, type(influence) targets(2) nlabel(0) name(plot_contract) replace nodraw
        local nr = r(n_rows)
        capture matrix PI = r(plot_data)
        local drc = _rc
        test_assert "[plots `family'] selected influence rows exposed" `drc'==0 & `nr'==60
        if `drc'==0 test_assert "[plots `family'] influence observation ids retained" PI[1,1]==2 & PI[1,2]==1
    }
    matrix PAFTER = e(b)
    test_assert "[plots `family'] preserves dataset and active estimates" _N==`original_n' & PAFTER[1,1]==PB[1,1]
    scalar plot_original_x = x1[1]
    replace x1 = x1+.01 in 1
    capture noisily `family'_dnn_plot, type(weights) view(kernel) name(plot_contract) replace nodraw
    test_assert "[plots `family'] snapshot kernel survives changed values" _rc==0
    if "`family'"!="np" {
        capture noisily `family'_dnn_plot, type(influence) name(plot_contract) replace nodraw
        test_assert "[plots `family'] snapshot influence survives changed values" _rc==0
    }
    replace x1 = scalar(plot_original_x) in 1
    // Restore exact original bytes before order check.
    estimates restore plot_`family'
    gsort -x1
    capture noisily `family'_dnn_plot, type(weights) view(kernel) name(plot_contract) replace nodraw
    test_assert "[plots `family'] snapshot kernel survives changed order" _rc==0
    capture noisily `family'_dnn_plot, targets(2) name(plot_contract) replace nodraw
    test_assert "[plots `family'] stored-result plots survive data changes" _rc==0
}

test_section "estimator weight plots: large-batch default target selection"
clear
set obs 12
gen double kw_x=(_n-6)/5
gen double kw_y=1+kw_x^2
gen double kw_den=.8+.05*kw_x
matrix KW_TARGETS=(-.8\-.6\-.4\-.2\0\.2\.4\.6\.8)
foreach family in cate np clate {
    local extra ""
    if "`family'"=="cate" local extra "pseudooutcome"
    if "`family'"=="clate" local extra "pseudooutcome denominator(kw_den)"
    quietly `family'_dnn kw_y, xcols(kw_x) x0matrix(KW_TARGETS) s(3) nose `extra'
    quietly `family'_dnn_plot, type(weights) view(kernel) name(plot_kernel_defaults) replace nodraw
    local nr=r(n_rows)
    local ids "`r(targets)'"
    test_assert "[plots `family'] large-batch kernel defaults to first/middle/last" `nr'==36 & "`ids'"=="1 5 9"
    quietly `family'_dnn_plot, type(weights) view(kernel) targets(8 2) name(plot_kernel_defaults) replace nodraw
    matrix KW_PLOT=r(plot_data)
    test_assert "[plots `family'] explicit kernel targets retain selection and order" rowsof(KW_PLOT)==24 & KW_PLOT[1,1]==8 & KW_PLOT[13,1]==2
    foreach view in profile cumulative {
        quietly `family'_dnn_plot, type(weights) view(`view') name(plot_kernel_defaults) replace nodraw
        local nr=r(n_rows)
        local ids "`r(targets)'"
        test_assert "[plots `family'] `view' defaults retain every target" `nr'==108 & "`ids'"=="1 2 3 4 5 6 7 8 9"
    }
}

test_section "estimator plot coherence: uncertainty and undefined target masks"
foreach family in cate np clate {
    clear
    set obs 20
    gen double x = _n
    gen double signal = 0
    gen double denominator = 1
    if "`family'"=="clate" replace signal = 2
    if "`family'"=="cate" quietly cate_dnn signal, xcols(x) x0(12) pseudooutcome
    if "`family'"=="np" quietly np_dnn signal, xcols(x) x0(12)
    if "`family'"=="clate" quietly clate_dnn signal, xcols(x) x0(12) denominator(denominator) pseudooutcome tuningmode(grid)
    quietly `family'_dnn_plot, name(plot_zero) replace nodraw
    matrix PZ = r(plot_data)
    if "`family'"=="clate" {
        test_assert "[plots CLATE] actual R fixture: zero SE retained with unavailable Wald bounds" PZ[1,6]==R_report_clate_zero_plot[1,1] & PZ[1,7]==R_report_clate_zero_plot[1,2] & missing(PZ[1,8]) & missing(PZ[1,9])
        quietly clate_dnn_plot, candidates(1) name(plot_zero_candidate) replace nodraw
        matrix PZC = r(plot_data)
        test_assert "[plots CLATE] candidate zero SE retained with unavailable Wald bounds" PZC[1,7]==0 & missing(PZC[1,8]) & missing(PZC[1,9])
        quietly clate_dnn_plot, type(tuning) name(plot_zero_tune) replace nodraw
        matrix PZT = r(plot_data)
        test_assert "[plots CLATE] tuning zero SE retained with unavailable Wald bounds" PZT[1,7]==0 & missing(PZT[1,8]) & missing(PZT[1,9])
    }
    else test_assert "[plots `family'] zero mean variance keeps degenerate interval" PZ[1,7]==0 & PZ[1,6]==0 & PZ[1,8]==0 & PZ[1,9]==0
    if "`family'"=="cate" quietly cate_dnn signal, xcols(x) x0(12) pseudooutcome nose
    if "`family'"=="np" quietly np_dnn signal, xcols(x) x0(12) nose
    if "`family'"=="clate" quietly clate_dnn signal, xcols(x) x0(12) denominator(denominator) pseudooutcome nose
    quietly `family'_dnn_plot, name(plot_nose) replace nodraw
    matrix PNO = r(plot_data)
    test_assert "[plots `family'] nose has missing uncertainty and bounds" missing(PNO[1,7]) & missing(PNO[1,8]) & missing(PNO[1,9])
    capture noisily `family'_dnn_plot, type(dashboard) style(bogus) replace nodraw
    test_assert "[plots `family'] dashboard rejects unknown style" _rc==198
    capture noisily `family'_dnn_plot, type(dashboard) style(profile) replace nodraw
    test_assert "[plots `family'] dashboard propagates invalid profile axis" _rc==198
    capture noisily `family'_dnn_plot, candidates(1) replace nodraw
    test_assert "[plots `family'] candidate request requires grid fit" _rc==498
    capture noisily `family'_dnn_plot, saving("build/coherence/unsupported.xyz") replace nodraw
    test_assert "[plots `family'] unsupported export extension rejected" _rc==198
}
clear
set obs 4
gen double x = _n
gen double sy = 1
gen double sd = cond(_n<=2,0,1)
matrix PPART = (0 \ 5)
quietly clate_dnn sy, denominator(sd) pseudooutcome xcols(x) x0matrix(PPART) s(4) nose
quietly clate_dnn_plot, style(forest) name(plot_partial) replace nodraw
local np = r(n_plotted)
local nr = r(n_rows)
local omitted = r(n_omitted)
matrix PP = r(plot_data)
test_assert "[plots CLATE] partial plot keeps legacy count and returns actual defined rows" `np'==2 & `nr'==1 & `omitted'==1 & PP[1,5]==2
replace sd = 0
quietly clate_dnn sy, denominator(sd) pseudooutcome xcols(x) x0matrix(PPART) s(4) nose
capture noisily clate_dnn_plot, name(plot_empty) replace nodraw
local rc = _rc
local nr = r(n_rows)
local np = r(n_plotted)
local omitted = r(n_omitted)
local matrices : r(matrices)
local hasdata : list posof "plot_data" in matrices
test_assert "[plots CLATE] empty plot remains usable without fabricated data rows" `rc'==0 & `nr'==0 & `np'==2 & `omitted'==2 & `hasdata'==0

// Instrument overlap uses the same selections with instrument-specific names.
rdnn_dgp_basic, n(50) seed(191)
gen byte a = runiform()<.5
gen byte d = a*(runiform()<.8)
gen double ep = .5
gen double md1 = .8
gen double my0 = x1
gen double my1 = x1+1.2
matrix PIV = (-1,0 \ 0,0 \ 1,0)
quietly clate_dnn y, xcols(x1 x2) x0matrix(PIV) dcol(d) acol(a) ecol(ep) my0col(my0) my1col(my1) md1col(md1) onesided
quietly clate_dnn_plot, type(overlap) targets(2 3) name(plot_iv) replace nodraw
local nr = r(n_rows)
matrix POVER = r(plot_data)
local names : colnames POVER
local armcol : list posof "instrument_1_mass" in names
test_assert "[plots CLATE] instrument overlap rows retain scientific labels and target IDs" `nr'==2 & `armcol'>0 & POVER[1,3]==2 & POVER[2,3]==3

test_section "estimator plot coherence: target-coordinate names cannot collide"
foreach coordinate in target_id treated_mass {
    rdnn_dgp_basic, n(40) seed(118)
    rename x1 `coordinate'
    matrix PCOLLISION = (-1,0 \ 0,0 \ 1,0)
    quietly cate_dnn y, xcols(`coordinate' x2) x0matrix(PCOLLISION) tcol(D) ecol(e) m0col(m0) m1col(m1)
    capture noisily cate_dnn_plot, type(overlap) targets(2 3) xindex(1) name(plot_collision) replace nodraw
    local rc = _rc
    matrix PC = r(plot_data)
    local names : colnames PC
    local coord : list posof "x0_1" in names
    test_assert "[plots coordinates] overlap handles covariate named `coordinate'" `rc'==0 & `coord'==1 & PC[1,1]==0 & PC[2,1]==1 & PC[1,3]==2
    capture noisily cate_dnn_plot, type(overlap) view(propensity) targets(2 3) name(plot_collision_auto) replace nodraw
    test_assert "[plots coordinates] overlap auto-axis handles `coordinate'" _rc==0
}

// Infeasible candidates retain missing status rows. Their missing values must
// not contaminate Wald bounds of available candidates (0 * missing is missing).
clear
set obs 20
gen double x = _n
gen double sy = sin(x)+x/10
gen double sd = 1
matrix P_INFEASIBLE = (5 \ 15)
quietly clate_dnn sy, denominator(sd) pseudooutcome xcols(x) ///
    x0matrix(P_INFEASIBLE) tuningmode(grid) sgrid(2 20) sdefault(2)
foreach kind in estimates tuning {
    quietly clate_dnn_plot, type(`kind') candidates(1) name(plot_infeasible) replace nodraw
    matrix P_AVAILABLE = r(plot_data)
    local cn : colnames P_AVAILABLE
    local sel : list posof "se" in cn
    local cil : list posof "ci_lower" in cn
    local ciu : list posof "ci_upper" in cn
    test_assert "[plots CLATE] `kind' valid candidate bounds survive infeasible grid rows" P_AVAILABLE[1,`sel']>0 & !missing(P_AVAILABLE[1,`cil'],P_AVAILABLE[1,`ciu']) & P_AVAILABLE[1,`ciu']>P_AVAILABLE[1,`cil']
}

// The finite n-s+k edge belongs to the exact rank table. Asymptotic
// weights can remain positive beyond that rank, including at n.
test_section "rank profile exact-support annotations respect the weight method"
clear
set obs 40
gen double support_plot_x=_n
gen double support_plot_y=sin(support_plot_x)
foreach method in exact asymp {
    quietly np_dnn support_plot_y, xcols(support_plot_x) x0(0) ///
        tuningmode(grid) sgrid(4 6) sdefault(4) method(`method') nose
    foreach selection in default candidate {
        local candidateopt ""
        if "`selection'"=="candidate" local candidateopt "candidates(1)"
        local svg "build/coherence/np-weights-support-`method'-`selection'.svg"
        quietly np_dnn_plot, type(weights) view(profile) `candidateopt' ///
            name(plot_support_contract) replace nodraw saving("`svg'")
        matrix SUPPORT_PROFILE=r(plot_data)
        if "`selection'"=="default" {
            local last_weight=SUPPORT_PROFILE[40,2]
            if "`method'"=="exact" test_assert "[plots exact] rank-table support ends at 37" `last_weight'==0
            else test_assert "[plots asymp] represented support extends beyond rank 37" `last_weight'>0
        }
        // Inspect the actual exported graph note, not a source-code pattern.
        tempname svgfile
        file open `svgfile' using "`svg'", read text
        local has_exact_edge 0
        file read `svgfile' line
        while r(eof)==0 {
            if strpos(`"`line'"',"Nominal exact-support edge: 37")>0 local has_exact_edge 1
            file read `svgfile' line
        }
        file close `svgfile'
        local expected=("`method'"=="exact" & "`selection'"=="default")
        test_assert "[plots `method'/`selection'] exported graph has only applicable exact-support note" `has_exact_edge'==`expected'
    }
}
