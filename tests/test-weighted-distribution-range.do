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
foreach engine in mata plugin {
    global rdistnn_distance_engine `engine'
    clear
    set obs 257
    gen double wd_x=cond(_n==1,129,cond(_n<=130,_n-2,_n-1))
    gen double wd_y=cond(mod(_n,2),1,2)
    replace wd_y=2^100 in 1
    quietly np_dnn wd_y, xcols(wd_x) x0(-1) k(1) s(128) nose
    matrix wd_diag=e(target_diagnostics)
    mata {
    wd_w=rdistnn_weights_exact(rdistnn_order(abs(st_data(.,"wd_x"):+1)),1,128)
    wd_z=st_data(.,"wd_y"):-1
    st_numscalar("wd_tail_weight",wd_w[1])
    st_numscalar("wd_sd_reference",sqrt(quadsum(wd_w:*wd_z:*wd_z)-quadsum(wd_w:*wd_z)^2))
    }
    test_assert "`engine': tiny positive tail remains in support" scalar(wd_tail_weight)>0
    test_assert "`engine': tiny-tail outcome preserves local SD" abs(wd_diag[1,colnumb(wd_diag,"weighted_sd_y")]-scalar(wd_sd_reference))<1e-15
    test_assert "`engine': requested distance engine ran" "$rdistnn_last_distance_engine"=="`engine'"
    gen long wd_order=cond(_n==1,258,_n)
    sort wd_order
    quietly np_dnn wd_y, xcols(wd_x) x0(-1) k(1) s(128) nose
    matrix wd_reordered=e(target_diagnostics)
    test_assert "`engine': outcome SD agrees after row permutation" abs(wd_reordered[1,colnumb(wd_reordered,"weighted_sd_y")]-scalar(wd_sd_reference))<1e-15
    gen double wd_neg=cond(mod(_n,2),1e-300,2e-300)
    replace wd_neg=(-5e307)*2 if wd_x==129
    quietly count if missing(wd_neg)
    local wd_neg_missing=r(N)
    test_assert "`engine': computed negative extreme is nonmissing" `wd_neg_missing'==0
    quietly np_dnn wd_neg, xcols(wd_x) x0(-1) k(1) s(128) nose
    matrix wd_negative=e(target_diagnostics)
    mata {
        wd_neg_y=st_data(.,"wd_neg")
        wd_neg_w=rdistnn_weights_exact(rdistnn_order(abs(st_data(.,"wd_x"):+1)),1,128)
        wd_neg_at=selectindex(st_data(.,"wd_x"):==129)
        st_numscalar("wd_neg_expected",-(wd_neg_y[wd_neg_at]*sqrt(wd_neg_w[wd_neg_at])))
    }
    test_assert "`engine': finite negative-extreme outcome SD" abs(wd_negative[1,colnumb(wd_negative,"weighted_sd_y")]/scalar(wd_neg_expected)-1)<2e-14
}
mata:
st_numscalar("wd_offset",rdistnn_weighted_sd((1e16\1e16+2\1e16\1e16+2),J(4,1,.25)))
st_numscalar("wd_tiny",rdistnn_weighted_sd((-1e-200\1e-200\1e300),(.5\.5\0))/1e-200)
wd_tail_w=(rdistnn_power2(-250)\.5\.5)
wd_tail_y=(rdistnn_power2(100)\1\2)
wd_tail_z=wd_tail_y:-1.5
wd_tail_reference=sqrt(quadsum(wd_tail_w:*wd_tail_z:*wd_tail_z)-quadsum(wd_tail_w:*wd_tail_z)^2)
st_numscalar("wd_tail",rdistnn_weighted_sd(wd_tail_y,wd_tail_w))
st_numscalar("wd_tail_expected",wd_tail_reference)
st_numscalar("wd_tail_reverse",rdistnn_weighted_sd(wd_tail_y[(3\2\1)],wd_tail_w[(3\2\1)]))
wd_tail_y[1]=-wd_tail_y[1]
st_numscalar("wd_tail_negative",rdistnn_weighted_sd(wd_tail_y,wd_tail_w))
st_numscalar("wd_q0",rdistnn_weighted_quantile((1::4),(0\.25\.5\.25),0))
st_numscalar("wd_q25",rdistnn_weighted_quantile((1::4),(0\.25\.5\.25),.25))
st_numscalar("wd_q75",rdistnn_weighted_quantile((1::4),(0\.25\.5\.25),.75))
st_numscalar("wd_q1",rdistnn_weighted_quantile((1::4),(0\.25\.5\.25),1))
end
test_assert "weighted SD preserves variation around a large offset" scalar(wd_offset)==1
test_assert "inactive large outcome does not erase tiny spread" abs(scalar(wd_tiny)-1)<2e-15
test_assert "positive tiny tail does not erase local spread" abs(scalar(wd_tail)-scalar(wd_tail_expected))<1e-15
test_assert "positive-tail SD agrees after reversing rows" abs(scalar(wd_tail_reverse)-scalar(wd_tail_expected))<1e-15
test_assert "negative tiny tail does not erase local spread" abs(scalar(wd_tail_negative)-scalar(wd_tail_expected))<1e-15
test_assert "zero probability preserves first-observation convention" scalar(wd_q0)==1
test_assert "quarter-mass boundary is exact" scalar(wd_q25)==2
test_assert "three-quarter boundary is exact" scalar(wd_q75)==3
test_assert "unit probability selects final positive mass" scalar(wd_q1)==4
macro drop rdistnn_distance_engine
do "tests/fixtures/weighted_cdf_fixtures.do"
