version 17
quietly do "src/_rdistnn_mata.ado"
do "tests/test_framework.do"
test_section "weighted diagnostic range and canonical CDF decisions"
foreach engine in mata plugin {
    global rdistnn_distance_engine `engine'
    clear
    set obs 8
    gen double wd_x=_n
    gen double wd_y=cond(mod(_n,2),-1e155,1e155)
    gen double wd_den=1
    quietly np_dnn wd_y, xcols(wd_x) x0(0) s(1) nose
    matrix wd_diag=e(target_diagnostics)
    test_assert "`engine': large finite outcome SD" abs(wd_diag[1,colnumb(wd_diag,"weighted_sd_y")]/1e155-1)<2e-15
    tempfile wd_saved
    estimates save `wd_saved'
    estimates clear
    estimates use `wd_saved'
    quietly rdistnn_predict, x0(0) post
    quietly np_dnn_diagnostics
    matrix wd_diag=r(diagnostics)
    test_assert "`engine': restored large finite outcome SD" abs(wd_diag[1,colnumb(wd_diag,"weighted_sd_y")]/1e155-1)<2e-15
    foreach magnitude in 1e200 1e-200 {
        replace wd_x=_n*`magnitude'
        replace wd_y=_n
        quietly cate_dnn wd_y, xcols(wd_x) x0(0) pseudooutcome s(2) nose
        matrix wd_diag=e(target_diagnostics)
        test_assert "`engine'/`magnitude': physical radius" abs(wd_diag[1,colnumb(wd_diag,"physical_radius_90")]/`magnitude'-6)<2e-14
        local wd_marker "$rdistnn_last_distance_engine"
        mata: wd_physical=rdistnn_physical_distances(st_data(.,"wd_x"),0)
        test_assert "physical distance helper preserves engine marker" "$rdistnn_last_distance_engine"=="`wd_marker'"
        quietly rdistnn_predict, x0(0) post
        quietly cate_dnn_diagnostics
        matrix wd_diag=r(diagnostics)
        test_assert "`engine'/`magnitude': rebuilt physical radius" abs(wd_diag[1,colnumb(wd_diag,"physical_radius_90")]/`magnitude'-6)<2e-14
        quietly clate_dnn wd_y, denominator(wd_den) pseudooutcome xcols(wd_x) x0(0) s(2) nose distance(l1)
        matrix wd_diag=e(target_diagnostics)
        test_assert "`engine'/`magnitude': CLATE physical radius" abs(wd_diag[1,colnumb(wd_diag,"physical_radius_90")]/`magnitude'-6)<2e-14
        quietly rdistnn_predict, x0(0) post
        matrix wd_diag=e(target_diagnostics)
        test_assert "`engine'/`magnitude': reused CLATE physical radius" abs(wd_diag[1,colnumb(wd_diag,"physical_radius_90")]/`magnitude'-6)<2e-14
    }
}
mata:
st_numscalar("wd_offset",rdistnn_weighted_sd((1e16\1e16+2\1e16\1e16+2),J(4,1,.25)))
st_numscalar("wd_tiny",rdistnn_weighted_sd((-1e-200\1e-200\1e300),(.5\.5\0))/1e-200)
st_numscalar("wd_q0",rdistnn_weighted_quantile((1::4),(0\.25\.5\.25),0))
st_numscalar("wd_q25",rdistnn_weighted_quantile((1::4),(0\.25\.5\.25),.25))
st_numscalar("wd_q75",rdistnn_weighted_quantile((1::4),(0\.25\.5\.25),.75))
st_numscalar("wd_q1",rdistnn_weighted_quantile((1::4),(0\.25\.5\.25),1))
end
test_assert "weighted SD preserves variation around a large offset" scalar(wd_offset)==1
test_assert "inactive large outcome does not erase tiny spread" abs(scalar(wd_tiny)-1)<2e-15
test_assert "zero probability preserves first-observation convention" scalar(wd_q0)==1
test_assert "quarter-mass boundary is exact" scalar(wd_q25)==2
test_assert "three-quarter boundary is exact" scalar(wd_q75)==3
test_assert "unit probability selects final positive mass" scalar(wd_q1)==4
macro drop rdistnn_distance_engine
do "tests/fixtures/weighted_cdf_fixtures.do"
