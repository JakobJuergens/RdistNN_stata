// Golden arrays regenerated from the live R reference; no R runtime needed.
do "tests/test_framework.do"
do "tests/fixtures/r_clate_fixtures.do"
clear
svmat double R_cl_input
rename (R_cl_input1 R_cl_input2 R_cl_input3 R_cl_input4 R_cl_input5 ///
    R_cl_input6 R_cl_input7 R_cl_input8 R_cl_input9 R_cl_input10) ///
    (x1 x2 Y D A e my0 my1 md0 md1)
test_section "clate_dnn: R reference estimates, covariance and diagnostics"
forvalues case=1/4 {
    local method "exact"
    local distance "l2"
    local k 1
    if inlist(`case',2,3) local k 2
    if `case'==3 local method "asymp"
    if `case'==4 local distance "mahalanobis"
    clate_dnn Y, xcols(x1 x2) x0matrix(R_cl_targets) ///
        dcol(D) acol(A) ecol(e) my0col(my0) my1col(my1) ///
        md0col(md0) md1col(md1) clip(.05) s(4) k(`k') ///
        method(`method') distance(`distance')
    matrix B = e(estimates)
    matrix V = e(theta_vcov)
    matrix C = e(components)
    matrix CV = e(component_vcov)
    matrix DG = e(diagnostics)
    mata: st_numscalar("cl_diff",max(abs(st_matrix("B")-st_matrix("R_cl_theta`case'"))))
    test_assert "R ratio case `case'" scalar(cl_diff)<1e-10
    mata: st_numscalar("cl_diff",max(abs(st_matrix("V")-st_matrix("R_cl_vcov`case'"))))
    test_assert "R paired cross-target covariance case `case'" scalar(cl_diff)<1e-8
    mata: st_numscalar("cl_diff",max(abs(st_matrix("C")-st_matrix("R_cl_components`case'"))))
    test_assert "R component means case `case'" scalar(cl_diff)<1e-10
    mata: st_numscalar("cl_diff",max(abs(st_matrix("e(se_targets)")-st_matrix("R_cl_se`case'"))))
    test_assert "R paired standard errors case `case'" scalar(cl_diff)<1e-9
    forvalues j=1/3 {
        local y = 2*`j'-1
        local d = 2*`j'
        test_assert "R YY covariance case `case' target `j'" abs(el(CV,`y',`y')-el(R_cl_paired`case',`j',1))<1e-9
        test_assert "R YD covariance case `case' target `j'" abs(el(CV,`y',`d')-el(R_cl_paired`case',`j',2))<1e-9
        test_assert "R DD covariance case `case' target `j'" abs(el(CV,`d',`d')-el(R_cl_paired`case',`j',3))<1e-9
    }
    mata: st_numscalar("cl_diff",max(abs(st_matrix("DG")[,(1..6,8..11,13,14)]-st_matrix("R_cl_diag`case'")[,1..12])))
    test_assert "R geometry and instrument diagnostics case `case'" scalar(cl_diff)<1e-10
    clate_dnn_diagnostics
    matrix CD = r(components)
    mata: st_numscalar("cl_diff",max(abs(st_matrix("CD")[,4..8]-st_matrix("R_cl_diag`case'")[,14..18])))
    test_assert "R first-stage uncertainty and paired component diagnostics case `case'" scalar(cl_diff)<1e-9
    // Test actual reference weight arrays, including stable ties.
    quietly weights_dnn x1 x2, x0matrix(R_cl_targets) k(`k') s(4) method(`method') distance(`distance')
    mata: st_matrix("CLW",st_matrix("r(weights)")[,1])
    mata: st_numscalar("cl_diff",max(abs(st_matrix("CLW")-st_matrix("R_cl_weights`case'")[,1])))
    test_assert "R weight array case `case'" scalar(cl_diff)<1e-12
}
clate_dnn Y, xcols(x1 x2) x0matrix(R_cl_targets) ///
    dcol(D) acol(A) ecol(e) my0col(my0) my1col(my1) ///
    md0col(md0) md1col(md1) clip(.05) tuningmode(grid) sgrid(3 5 3) sdefault(3)
matrix GB = e(tuning_estimates)
matrix GS = e(tuning_se)
forvalues j=1/3 {
    test_assert "R grid s3 target `j'" abs(el(GB,1,`j')-el(R_cl_grid3,1,`j'))<1e-10
    test_assert "R grid SE s3 target `j'" abs(el(GS,1,`j')-el(R_cl_grid3,2,`j'))<1e-9
    test_assert "R grid s5 target `j'" abs(el(GB,2,`j')-el(R_cl_grid5,1,`j'))<1e-10
    test_assert "R grid SE s5 target `j'" abs(el(GS,2,`j')-el(R_cl_grid5,2,`j'))<1e-9
}
