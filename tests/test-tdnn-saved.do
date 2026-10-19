// Two-scale fitted results must survive data changes and numeric serialization.
version 17
do "tests/test_framework.do"
test_section "TDNN saved reporting and prediction"
clear
set obs 24
gen double ts_x1=(_n-12)/5
gen double ts_x2=sin(_n)
gen double ts_y=2+ts_x1^2-.3*ts_x2
matrix TS_A=(-.5,0 \ .5,.2)
matrix TS_B=(-.2,.1 \ .8,-.1 \ 0,0)
foreach method in exact asymp {
    quietly np_tdnn ts_y, xcols(ts_x1 ts_x2) x0matrix(TS_A) ///
        s1(3) s2(7) method(`method') distance(l2_scaled) scalerule(iqr)
    test_assert "`method' initial covariance method is signed" "`e(se_method)'"=="signed_rank_shift"
    estimates store ts_original
    matrix TS_fit=e(estimates)
    matrix TS_cov=e(estimate_vcov)
    quietly rdistnn_results, type(vcov)
    matrix TS_out=r(vcov)
    mata: st_numscalar("ts_error",max(abs(st_matrix("TS_out")-st_matrix("TS_cov"))))
    test_assert "`method' reporting returns full covariance" scalar(ts_error)<1e-12
    quietly rdistnn_results, type(weights)
    matrix TS_weights=r(data)
    test_assert "`method' weight export has every observation and target" rowsof(TS_weights)==48 & colsof(TS_weights)==5
    quietly rdistnn_predict, x0matrix(TS_B) sefit store(ts_derived)
    matrix TS_prediction=r(fit)
    matrix TS_predictionV=r(vcov)
    mata: st_numscalar("ts_error",max(abs(st_matrix("e(estimates)")-st_matrix("TS_fit"))))
    test_assert "`method' prediction preserves active estimates" scalar(ts_error)==0 & e(m)==2
    quietly np_tdnn ts_y, xcols(ts_x1 ts_x2) x0matrix(TS_B) ///
        s1(3) s2(7) method(`method') distance(l2_scaled) scalerule(iqr)
    mata: st_numscalar("ts_error",max(abs(vec(st_matrix("TS_prediction"))-vec(st_matrix("e(estimates)")))))
    test_assert "`method' prediction equals fresh batch fit" scalar(ts_error)<1e-12
    mata: st_numscalar("ts_error",max(abs(st_matrix("TS_predictionV")-st_matrix("e(estimate_vcov)"))))
    test_assert "`method' predicted covariance equals fresh batch fit" scalar(ts_error)<1e-12
    estimates restore ts_original
    tempfile ts_saved
    quietly estimates save "`ts_saved'", replace
    preserve
    clear
    quietly estimates use "`ts_saved'"
    quietly rdistnn_predict, x0matrix(TS_B) sefit post
    mata: st_numscalar("ts_error",max(abs(st_matrix("e(estimate_vcov)")-st_matrix("TS_predictionV"))))
    test_assert "`method' cold saved model predicts without active data" scalar(ts_error)<1e-12 & e(m)==3
    test_assert "`method' posted covariance method is signed" "`e(se_method)'"=="signed_rank_shift"
    quietly np_tdnn_diagnostics
    test_assert "`method' derived fit has signed diagnostics" r(n_rows)==3
    quietly np_tdnn_plot, type(estimates) nodraw replace
    test_assert "`method' derived fit produces estimate plot data" rowsof(r(plot_data))==3
    quietly np_tdnn_plot, type(weights) nodraw replace
    capture rdistnn_candidate, candidate(1)
    test_assert "TDNN reports absent candidate grid" _rc==498
    restore
    estimates restore ts_derived
    test_assert "stored derived fit retains signed method" "`e(se_method)'"=="signed_rank_shift"
    test_assert "stored derived fit retains two-scale controls" e(s1)==3 & e(s2)==7 & "`e(cmd)'"=="np_tdnn"
}
quietly np_tdnn ts_y, xcols(ts_x1 ts_x2) x0matrix(TS_A) s1(3) s2(24) nose
quietly rdistnn_predict, x0matrix(TS_B)
test_assert "nose predictions leave uncertainty missing" missing(el(r(se_fit),1,1))
quietly rdistnn_predict, x0matrix(TS_B) post
test_assert "posted nose TDNN has no covariance method" "`e(se_method)'"=="none"
capture rdistnn_predict, x0matrix(TS_B) sefit
test_assert "saved boundary fit rejects infeasible covariance" _rc==198
