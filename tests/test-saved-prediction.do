// Saved prediction and serialization; reference arithmetic remains local.
version 17
do "tests/test_framework.do"
test_section "saved model snapshots for every estimator"
clear
set obs 12
gen double fr_x=_n/4
gen double fr_y=sin(fr_x)+fr_x
gen double fr_den=1+fr_x/10
gen double fr_d=2*fr_y
gen double fr_zero=0
gen double fr_positive=1
foreach family in np cate clate plqr {
    if "`family'"=="np" quietly np_dnn fr_y, xcols(fr_x) x0(0) s(3)
    if "`family'"=="cate" quietly cate_dnn fr_y, pseudooutcome xcols(fr_x) x0(0) s(3)
    if "`family'"=="clate" quietly clate_dnn fr_y, pseudooutcome denominator(fr_den) xcols(fr_x) x0(0) s(3)
    if "`family'"=="plqr" quietly plqr_dnn fr_positive, x(fr_x) x0(0) d(fr_d) md(fr_zero) g(fr_zero) thetaeval(0) s(3)
    test_assert "`family': snapshot format is versioned" e(model_version)==4
    test_assert "`family': untied localized fit records zero notice counts" e(n_tied_targets)==0 & e(n_weak_targets)==0 & r(n_tied_targets)==0 & r(n_weak_targets)==0
    matrix FR_NOTICES=e(target_notices)
    foreach field in inputs X signals ranks distances groups row_ids {
        local model_matrices : e(matrices)
        local required "model_`field'"
        local snapshot_present : list required in model_matrices
        capture matrix FR_SNAPSHOT=e(model_`field')
        test_assert "`family': e(model_`field') is stored" `snapshot_present'==1
        if `snapshot_present' {
            test_assert "`family': e(model_`field') retains the full fitted sample" rowsof(FR_SNAPSHOT)==e(N)
        }
    }
    if e(model_version)==4 {
        matrix FR_X=e(model_X)
        matrix FR_SIGNALS=e(model_signals)
        tempfile serialized
        quietly estimates save `serialized', replace
        preserve
        clear
        quietly estimates use `serialized'
        mata: st_numscalar("_fr_cold",max((max(abs(st_matrix("FR_X")-st_matrix("e(model_X)"))),max(abs(st_matrix("FR_SIGNALS")-st_matrix("e(model_signals)"))))))
        test_assert "`family': .ster restore retains numeric model state without data" scalar(_fr_cold)==0
        mata: st_numscalar("_fr_notices",all(st_matrix("FR_NOTICES"):==st_matrix("e(target_notices)")))
        test_assert "`family': .ster restore retains notice facts without data" scalar(_fr_notices) & e(n_tied_targets)==0 & e(n_weak_targets)==0
        restore
    }
}

test_section "PLQR near-tie geometry agrees with the R distance convention"
do tests/fixtures/r_reuse_fixtures.do
preserve
clear
set obs 3
gen double tie_x1=cond(_n<3,1,2)
gen double tie_x2=cond(_n==1,2^-26,0)
gen double tie_y=10
gen double tie_d=cond(_n==2,3,1)
gen double tie_zero=0
quietly plqr_dnn tie_y, x(tie_x1 tie_x2) x0(0 1) d(tie_d) md(tie_zero) g(tie_zero) thetaeval(0) s(2)
quietly rdistnn_predict, x0(0 0) sefit
matrix FR_TIE=r(fit)
quietly plqr_dnn tie_y, x(tie_x1 tie_x2) x0(0 0) d(tie_d) md(tie_zero) g(tie_zero) thetaeval(0) s(2)
matrix FR_TIEFRESH=e(moment)
// Equal rounded distances average total mass1 across the two represented
// nearest points. Their scores are3/2 and1/2, giving independent moment1.
test_assert "PLQR fresh fit and saved prediction resolve near ties identically to R" abs(FR_TIE[1,1]-1)<1e-12 & abs(FR_TIEFRESH[1,1]-R_reuse_plqr_neartie[1,1])<1e-12
restore

test_section "all estimators: cold predictions, engines, exports and failures"
matrix FR_NEW=(.25\1.25\2.75)
gen double fr_d2=cos(fr_x)
foreach method in exact asymp {
    foreach family in np cate clate plqr {
        local opts "xcols(fr_x) x0matrix(FR_NEW) s(4) k(2) method(`method')"
        if "`family'"=="np" quietly np_dnn fr_y, `opts'
        if "`family'"=="cate" quietly cate_dnn fr_y, pseudooutcome `opts'
        if "`family'"=="clate" quietly clate_dnn fr_y, pseudooutcome denominator(fr_den) `opts'
        if "`family'"=="plqr" quietly plqr_dnn fr_y, x(fr_x) x0(.25) d(fr_d fr_d2) md(fr_zero fr_zero) g(fr_zero) thetaeval(.3 -.2) s(4) k(2) method(`method')
        tempname saved
        estimates store `saved'
        matrix FR_REF=e(estimates)
        matrix FR_REFV=e(estimate_vcov)
        matrix FR_RANKS=e(model_ranks)
        foreach engine in mata naive plugin {
            global rdistnn_jackknife_engine `engine'
            local targetopt "x0matrix(FR_NEW)"
            if "`family'"=="plqr" local targetopt "x0(.25)"
            quietly rdistnn_predict, `targetopt' sefit
            test_assert "`family' `method' `engine': derived batch returns its retained tie count" r(n_tied_targets)==e(n_tied_targets) & r(n_weak_targets)==e(n_weak_targets)
            matrix FR_GOT=r(fit)
            matrix FR_GOTV=r(vcov)
            mata: st_numscalar("_fr_error",max((max(abs(vec(st_matrix("FR_GOT"))'-st_matrix("FR_REF"))),max(abs(st_matrix("FR_GOTV")-st_matrix("FR_REFV"))))))
            test_assert "`family' `method' `engine': saved evaluation matches fresh estimate and covariance" scalar(_fr_error)<1e-10
        }
        macro drop rdistnn_jackknife_engine
        tempfile cold tablefile
        quietly estimates save `cold', replace
        preserve
        clear
        mata: mata clear
        quietly estimates use `cold'
        quietly rdistnn_predict, `targetopt' sefit store(fr_derived)
        matrix FR_COLD=r(fit)
        mata: st_numscalar("_fr_error",max(abs(vec(st_matrix("FR_COLD"))'-st_matrix("FR_REF"))))
        test_assert "`family' `method': cold .ster prediction uses no Mata globals or active data" scalar(_fr_error)<1e-10 & _N==0
        quietly rdistnn_results, saving(`tablefile') replace
        matrix FR_REPORT=r(data)
        matrix FR_INDEX=r(index)
        local q=cond("`family'"=="plqr",2,1)
        test_assert "`family' `method': graph-free report retains target/component rows" rowsof(FR_REPORT)==e(m)*`q' & colsof(FR_INDEX)==2
        quietly rdistnn_results, type(vcov) saving(`tablefile') replace
        matrix FR_COVDATA=r(data)
        test_assert "`family' `method': covariance export includes explicit component and target IDs" rowsof(FR_COVDATA)==(e(m)*`q')^2 & colsof(FR_COVDATA)==5
        quietly rdistnn_results, type(weights)
        matrix FR_W=r(data)
        test_assert "`family' `method': saved weight table keeps every fitted observation" rowsof(FR_W)==e(N)*e(m)
        matrix FR_BAD=(1,2)
        capture rdistnn_predict, x0matrix(FR_BAD) post
        local badrc=_rc
        matrix FR_AFTER=e(estimates)
        mata: st_numscalar("_fr_error",max(abs(st_matrix("FR_AFTER")-st_matrix("FR_REF"))))
        test_assert "`family' `method': prediction error preserves caller and active fit" `badrc'==198 & _N==0 & scalar(_fr_error)==0
        quietly estimates restore fr_derived
        test_assert "`family' `method': stored derived fit retains its numeric model" e(model_version)==4
        estimates drop fr_derived
        restore
        quietly estimates restore `saved'
        estimates drop `saved'
    }
}

test_section "PLQR batch moments retain supplied coefficients separately"
quietly plqr_dnn fr_y, x(fr_x) x0(.25) d(fr_d fr_d2) md(fr_zero fr_zero) g(fr_zero) thetaeval(.3 -.2) s(4) k(2)
tempname plparent
estimates store `plparent'
quietly rdistnn_predict, x0matrix(FR_NEW) sefit post
matrix FR_PLQR=r(fit)
matrix FR_THETA=r(theta_eval)
test_assert "PLQR predictions are component-by-target moments" rowsof(FR_PLQR)==2 & colsof(FR_PLQR)==3 & "`r(quantity)'"=="moment"
test_assert "PLQR supplied theta is retained separately" FR_THETA[1,1]==.3 & FR_THETA[1,2]==-.2 & "`e(inference_target)'"=="moment"
quietly rdistnn_results, targets(3 1)
matrix FR_SELECTED=r(fit)
matrix FR_INDEX=r(index)
test_assert "PLQR selection preserves original target and component IDs" FR_INDEX[1,1]==3 & FR_INDEX[2,2]==2 & FR_INDEX[3,1]==1
forvalues t=1/3 {
    quietly plqr_dnn fr_y, x(fr_x) x0(`=FR_NEW[`t',1]') d(fr_d fr_d2) md(fr_zero fr_zero) g(fr_zero) thetaeval(.3 -.2) s(4) k(2)
    matrix FR_SINGLE=e(moment)
    test_assert "PLQR batch target `t' matches original single-target evaluator" abs(FR_SINGLE[1,1]-FR_PLQR[1,`t'])<1e-12 & abs(FR_SINGLE[1,2]-FR_PLQR[2,`t'])<1e-12
}
quietly estimates restore `plparent'
estimates drop `plparent'

test_section "undefined CLATE and unavailable covariance keep full identities"
gen double fr_partial=(_n>1)
matrix FR_ENDPOINTS=(0\3)
quietly clate_dnn fr_y, denominator(fr_partial) pseudooutcome xcols(fr_x) x0matrix(FR_ENDPOINTS) s(12) nose
quietly rdistnn_predict, x0matrix(FR_ENDPOINTS) post
matrix FR_RATIO=r(fit)
matrix FR_SE=r(se_fit)
test_assert "CLATE prediction keeps undefined ratios and nose uncertainty" missing(FR_RATIO[1,1]) & FR_RATIO[1,2]<. & missing(FR_SE[1,2])
quietly rdistnn_results, type(status) targets(2 1)
matrix FR_STATUS=r(data)
test_assert "CLATE status selection retains original IDs and undefined status" FR_STATUS[1,1]==2 & FR_STATUS[2,1]==1 & FR_STATUS[2,2]==1
quietly rdistnn_results, type(components)
matrix FR_COMPONENTS=r(data)
test_assert "CLATE undefined ratio retains its zero first-stage component" FR_COMPONENTS[1,3]==0
capture rdistnn_predict, x0(0) sefit post
test_assert "unavailable deletion covariance is rejected without changing nose fit" _rc==198 & e(s)==12 & "`e(se_method)'"=="none"

test_section "R-reference prediction and full target covariance"
do tests/fixtures/r_reuse_fixtures.do
foreach method in exact asymp {
    foreach family in np cate clate plqr {
        local opts "xcols(fr_x) x0(0) s(4) k(2) method(`method')"
        if "`family'"=="np" quietly np_dnn fr_y, `opts'
        if "`family'"=="cate" quietly cate_dnn fr_y, pseudooutcome `opts'
        if "`family'"=="clate" quietly clate_dnn fr_y, pseudooutcome denominator(fr_den) `opts'
        if "`family'"=="plqr" quietly plqr_dnn fr_y, x(fr_x) x0(0) d(fr_d fr_d2) md(fr_zero fr_zero) g(fr_zero) thetaeval(.3 -.2) s(4) k(2) method(`method')
        quietly rdistnn_predict, x0matrix(FR_NEW) sefit
        matrix FR_FIT=r(fit)
        matrix FR_V=r(vcov)
        mata: st_numscalar("_fr_error",max((max(abs(st_matrix("FR_FIT")-st_matrix("R_reuse_`family'_`method'_fit"))),max(abs(st_matrix("FR_V")-st_matrix("R_reuse_`family'_`method'_vcov"))))))
        test_assert "`family' `method': prediction and joint covariance match independent R frozen deletions" scalar(_fr_error)<1e-10
    }
}

test_section "fixed singular Mahalanobis geometry and error preservation"
gen double fr_x2=2*fr_x
matrix FR_MNEW=(.25,.5\1.25,2.5\2.75,5.5)
quietly np_dnn fr_y, xcols(fr_x fr_x2) x0matrix(FR_MNEW) s(4) distance(mahalanobis)
matrix FR_REF=e(estimates)
matrix FR_REFV=e(estimate_vcov)
quietly np_dnn fr_y, xcols(fr_x fr_x2) x0(0 0) s(4) distance(mahalanobis)
preserve
replace fr_y=1e6
gsort -fr_x
quietly rdistnn_predict, x0matrix(FR_MNEW) sefit
matrix FR_FIT=r(fit)
matrix FR_V=r(vcov)
mata: st_numscalar("_fr_error",max((max(abs(st_matrix("FR_FIT")-st_matrix("FR_REF"))),max(abs(st_matrix("FR_V")-st_matrix("FR_REFV"))))))
test_assert "prediction uses fixed fitted geometry and signals after active inputs change" scalar(_fr_error)<1e-10 & fr_y[1]==1e6
capture rdistnn_results, saving("build/nonexistent-fr-folder/report.dta")
test_assert "failed report save preserves data and fit" _rc!=0 & _N==12 & fr_y[1]==1e6 & e(s)==4
restore
