// One inferential/reporting contract for scalar CATE, CLATE and NP fits.
version 17
do "tests/test_framework.do"
do "tests/fixtures/r_reporting_fixtures.do"
clear
set obs 24
gen double x = (_n-12)/4
gen double signal = 2+sin(x)+cos(_n)/5
gen double unit = 1
matrix REPORT_X0 = (-1 \ 0 \ 1)
test_section "estimator coherence: actual summary tables and replay"
foreach family in cate clate np {
    local opts ""
    if "`family'"=="cate" local opts "pseudooutcome"
    if "`family'"=="clate" local opts "pseudooutcome denominator(unit)"
    quietly `family'_dnn signal, xcols(x) x0matrix(REPORT_X0) s(4) `opts'
    matrix RT = r(table)
    matrix RS = r(summary)
    if rowsof(RS)==3 & colsof(RS)==6 {
        mata: st_numscalar("report_r_error",max(abs(st_matrix("RS")-st_matrix("R_report_`family'95"))))
        test_assert "`family' summary agrees with live R reference" scalar(report_r_error)<1e-10
    }
    matrix RB = e(b)
    matrix RV = e(V)
    local rnames : rownames RT
    local snames : colnames RS
    test_assert "`family' standard table is actually 9 by m" rowsof(RT)==9 & colsof(RT)==3
    test_assert "`family' common summary has six columns" rowsof(RS)==3 & colsof(RS)==6
    test_assert "`family' standard row names" "`rnames'"=="b se z pvalue ll ul df crit eform"
    test_assert "`family' common column names" "`snames'"=="estimate std_error z pvalue ci_lower ci_upper"
    if rowsof(RT)==9 & rowsof(RS)==3 & colsof(RS)==6 {
        mata: st_numscalar("report_error",max(abs(st_matrix("RT")[1..6,.]'-st_matrix("RS"))))
        test_assert "`family' standard and tidy statistics agree" scalar(report_error)<1e-12
        test_assert "`family' estimates and standard errors follow fit" abs(RT[1,2]-RB[1,2])<1e-12 & abs(RT[2,2]-sqrt(RV[2,2]))<1e-12
        test_assert "`family' z and p use normal inference" abs(RT[3,2]-RB[1,2]/sqrt(RV[2,2]))<1e-12 & abs(RT[4,2]-2*normal(-abs(RT[3,2])))<1e-12
    }
    matrix CE = e(estimates)
    matrix CS = e(se_targets)
    matrix CV = e(estimate_vcov)
    test_assert "`family' canonical estimates/SEs preserve all targets" colsof(CE)==3 & colsof(CS)==3 & !missing(CS[1,2])
    test_assert "`family' canonical covariance preserves batch dimensions" rowsof(CV)==3 & colsof(CV)==3 & abs(CV[1,2]-RV[1,2])<1e-12
    test_assert "`family' ordered input signature stored" "`e(ordered_signature)'"!="" & "`e(signature_vars)'"!=""
    return clear
    quietly `family'_dnn, level(90)
    matrix R90 = r(table)
    matrix S90 = r(summary)
    local replay_level = r(level)
    if rowsof(S90)==3 & colsof(S90)==6 {
        mata: st_numscalar("report_r_error",max(abs(st_matrix("S90")-st_matrix("R_report_`family'90"))))
        test_assert "`family' replay agrees with R at 90 percent" scalar(report_r_error)<1e-10
    }
    test_assert "`family' replay rebuilds actual standard table" rowsof(R90)==9 & colsof(R90)==3
    test_assert "`family' replay records level" `replay_level'==90
    if rowsof(R90)==9 {
        test_assert "`family' replay level changes confidence limits" abs(R90[5,2]-(RB[1,2]-invnormal(.95)*sqrt(RV[2,2])))<1e-12 & abs(R90[8,2]-invnormal(.95))<1e-12
    }
    mata: st_numscalar("report_error",max(abs(st_matrix("e(b)")-st_matrix("RB")))+max(abs(st_matrix("e(V)")-st_matrix("RV"))))
    test_assert "`family' replay does not mutate fit" scalar(report_error)==0
    estimates store report_saved
    quietly np_dnn signal, xcols(x) x0(0) s(3) nose
    estimates restore report_saved
    return clear
    quietly `family'_dnn
    matrix RRESTORE = r(table)
    test_assert "`family' restored summary remains full target" rowsof(RRESTORE)==9 & colsof(RRESTORE)==3

    quietly `family'_dnn signal, xcols(x) x0matrix(REPORT_X0) s(4) `opts' nose
    matrix RN = r(table)
    matrix SN = r(summary)
    test_assert "`family' nose preserves estimates without precision" rowsof(RN)==9 & !missing(RN[1,2]) & missing(RN[2,2]) & missing(RN[5,2])
    matrix CN = e(se_targets)
    test_assert "`family' nose canonical SE is missing" colsof(CN)==3 & missing(CN[1,2])
}

test_section "estimator coherence: CLATE exceptional targets"
clear
set obs 4
gen double x=_n
gen double signal=cond(_n==1,2,cond(_n==2,1,cond(_n==3,4,3)))
gen double stage=cond(_n<=2,0,cond(_n==3,-1,2))
matrix REPORT_PARTIAL=(0 \ 5)
quietly clate_dnn signal, denominator(stage) pseudooutcome xcols(x) x0matrix(REPORT_PARTIAL) s(3)
matrix RP = r(table)
matrix SP = r(summary)
test_assert "partial CLATE standard table retains undefined target" rowsof(RP)==9 & colsof(RP)==2 & missing(RP[1,1]) & abs(RP[1,2]-2.6)<1e-12
test_assert "partial CLATE summary preserves target order" rowsof(SP)==2 & colsof(SP)==6 & missing(SP[1,1]) & abs(SP[2,1]-2.6)<1e-12
replace stage=0
quietly clate_dnn signal, denominator(stage) pseudooutcome xcols(x) x0matrix(REPORT_PARTIAL) s(3)
matrix RZ = r(table)
test_assert "all undefined CLATE reports without e(b)" rowsof(RZ)==9 & colsof(RZ)==2 & missing(RZ[1,1]) & missing(RZ[1,2])
replace stage=1
replace signal=2
quietly clate_dnn signal, denominator(stage) pseudooutcome xcols(x) x0(0) s(2)
matrix RZERO = r(table)
matrix RZS = r(summary)
mata: st_numscalar("report_zero_error", max(abs(st_matrix("RZS")[1,1..2]-st_matrix("R_report_clate_zero_summary")[1,1..2])))
test_assert "CLATE zero SE summary matches actual R status convention" scalar(report_zero_error)<1e-12 & missing(RZS[1,3]) & missing(RZS[1,6])
test_assert "CLATE zero SE retains value but masks Wald fields" rowsof(RZERO)==9 & RZERO[1,1]==2 & RZERO[2,1]==0 & missing(RZERO[3,1]) & missing(RZERO[4,1]) & missing(RZERO[5,1])

// Ordinary mean/CATE zero variance follows its R reporting convention rather
// than applying the affine-ratio status mask to another estimator.
foreach family in cate np {
    local opts ""
    if "`family'"=="cate" local opts "pseudooutcome"
    quietly `family'_dnn signal, xcols(x) x0(0) s(2) `opts'
    matrix RMZERO = r(table)
    test_assert "`family' zero-variance mean bounds retained" rowsof(RMZERO)==9 & RMZERO[1,1]==2 & RMZERO[2,1]==0 & RMZERO[5,1]==2 & RMZERO[6,1]==2
}
