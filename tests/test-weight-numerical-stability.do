// Stable tail masses and adjacent-scale TDNN coefficients through shared tables.
do "tests/test_framework.do"
quietly run "src/_rdistnn_mata.ado"
quietly run "src/_rdistnn_tdnn_mata.ado"
test_section "Weight stability: analytic Beta(1,b) tiny tails"
mata {
    WNS_w=rdistnn_weights_asymp((1::100),1,80)
    WNS_err=0
    for(WNS_r=50;WNS_r<=100;WNS_r++) {
        WNS_lo=(WNS_r-1)/100; WNS_hi=WNS_r/100
        WNS_ref=exp(80*ln1p(-WNS_lo))
        if(WNS_hi<1) WNS_ref=WNS_ref*(-expm1(80*(ln1p(-WNS_hi)-ln1p(-WNS_lo))))
        WNS_err=max((WNS_err,abs(WNS_w[WNS_r]/WNS_ref-1)))
    }
    st_numscalar("WNS_error",WNS_err)
    st_numscalar("WNS_positive",all(WNS_w:>0))
}
test_assert "representable tiny masses retain exact support" scalar(WNS_positive)==1
test_assert "upper-tail masses match analytic powers relatively" scalar(WNS_error)<2e-11

test_section "Weight stability: large bias dimension with adjacent scales"
mata {
    WNS_td=rdistnn_tdnn_table(100,79,80,1e17,"asymp",0)
    WNS_c=-1/expm1((2/1e17)*ln1p(1/79))
    WNS_a=rdistnn_weights_asymp((1::100),1,79)
    WNS_b=rdistnn_weights_asymp((1::100),1,80)
    WNS_ref=WNS_c*WNS_a+(1-WNS_c)*WNS_b
    st_numscalar("WNS_finite",!hasmissing(WNS_td))
    st_numscalar("WNS_error",max(abs(WNS_td-WNS_ref))/max(abs(WNS_ref)))
}
test_assert "adjacent-scale TDNN remains finite at large dim" scalar(WNS_finite)==1
test_assert "adjacent-scale TDNN matches stable coefficient oracle" scalar(WNS_error)<2e-12

test_section "Weight stability: reflected tails, endpoints and ordinary shapes"
mata {
    WNS_ref=1.6542999783940550651810050759970065291673e-68
    st_numscalar("WNS_error",max((abs(rdistnn_beta_interval(1,100,.79,.8)/WNS_ref-1),
        abs(rdistnn_beta_interval(100,1,.2,.21)/WNS_ref-1))))
    WNS_error2=0
    for(WNS_a=1;WNS_a<=7;WNS_a++) {
        WNS_error2=max((WNS_error2,abs(rdistnn_beta_interval(WNS_a,8-WNS_a,0,1)-1)))
        WNS_error2=max((WNS_error2,abs(rdistnn_beta_interval(WNS_a,8-WNS_a,0,0))))
        WNS_error2=max((WNS_error2,abs(rdistnn_beta_interval(WNS_a,8-WNS_a,1,1))))
    }
    st_numscalar("WNS_endpoints",WNS_error2)
    WNS_coef=rdistnn_tdnn_coefficients(1e9,1e9+1,3)
    st_numscalar("WNS_coef_error",abs(WNS_coef[1]/(-1500000000.25)-1))
    WNS_error2=0
    for(WNS_k=1;WNS_k<=7;WNS_k++) {
        WNS_w=rdistnn_weights_asymp((24::1),WNS_k,7)
        for(WNS_i=1;WNS_i<=24;WNS_i++) {
            WNS_r=25-WNS_i; WNS_ref=0
            for(WNS_a=1;WNS_a<=WNS_k;WNS_a++) {
                WNS_ref=WNS_ref+(ibeta(WNS_a,8-WNS_a,WNS_r/24)-ibeta(WNS_a,8-WNS_a,(WNS_r-1)/24))/WNS_k
            }
            WNS_error2=max((WNS_error2,abs(WNS_w[WNS_i]-WNS_ref)))
        }
        WNS_error2=max((WNS_error2,abs(sum(WNS_w)-1)))
    }
    st_numscalar("WNS_ordinary",WNS_error2)
}
test_assert "direct and reflected tails match independent decimal fixture" scalar(WNS_error)<2e-12
test_assert "full and zero-width endpoint intervals are exact" scalar(WNS_endpoints)==0
test_assert "adjacent scales match independent decimal coefficient" scalar(WNS_coef_error)<3e-16
test_assert "all k and reversed ranks preserve ordinary weights and unit mass" scalar(WNS_ordinary)<2e-14

test_section "Weight stability: public scalar and batched weight consumers"
clear
set obs 100
gen double wns_x=_n-1
gen double wns_y=sin(_n/10)
gen double wns_den=1+_n/200
matrix WNS_X0=(0\0\2)
quietly weights_dnn wns_x, x0matrix(WNS_X0) k(1) s(80) method(asymp)
matrix WNS_W=r(weights)
mata: st_numscalar("WNS_positive",all(st_matrix("WNS_W"):>0))
test_assert "public batch DNN weights retain every positive tail" scalar(WNS_positive)==1
quietly weights_tdnn wns_x, x0matrix(WNS_X0) s1(79) s2(80) dim(100000000000000000) method(asymp) nonormalize
matrix WNS_T=r(rank_weights)
mata: st_numscalar("WNS_finite",!hasmissing(st_matrix("WNS_T")))
test_assert "public TDNN nonormalize retains finite large-dim weights" scalar(WNS_finite)==1
quietly np_dnn wns_y, xcols(wns_x) x0matrix(WNS_X0) s(80) method(asymp)
matrix WNS_np=e(estimates)
mata: st_numscalar("WNS_mean",rdistnn_weights_asymp((1::100),1,80)'*st_data(.,"wns_y"))
scalar WNS_error=max(abs(el(WNS_np,1,1)-WNS_mean),abs(el(WNS_np,1,2)-WNS_mean))
test_assert "NP duplicate targets consume the shared stable table" scalar(WNS_error)<2e-12
quietly rdistnn_predict, x0matrix(WNS_X0) sefit
matrix WNS_prediction=r(fit)
mata: st_numscalar("WNS_error",max(abs(st_matrix("WNS_prediction")-st_matrix("WNS_np"))))
test_assert "saved NP prediction retains stable asymptotic weights" scalar(WNS_error)<2e-12
quietly cate_dnn wns_y, xcols(wns_x) x0(0) pseudooutcome s(80) method(asymp)
scalar WNS_cate=el(e(b),1,1)
test_assert "CATE pseudooutcome uses the same tail-preserving weights" abs(scalar(WNS_cate)-el(WNS_np,1,1))<2e-12
quietly clate_dnn wns_y, xcols(wns_x) x0(0) pseudooutcome denominator(wns_den) s(80) method(asymp)
scalar WNS_clate=el(e(b),1,1)
mata: st_numscalar("WNS_ratio",(rdistnn_weights_asymp((1::100),1,80)'*st_data(.,"wns_y"))/(rdistnn_weights_asymp((1::100),1,80)'*st_data(.,"wns_den")))
test_assert "CLATE ratio uses the same stable component weights" abs(scalar(WNS_clate)-scalar(WNS_ratio))<2e-12
foreach normopt in "" "nonormalize" {
    quietly weights_tdnn wns_x, x0matrix(WNS_X0) s1(79) s2(80) method(asymp) `normopt'
    matrix WNS_TW=r(weights)
    quietly np_tdnn wns_y, xcols(wns_x) x0matrix(WNS_X0) s1(79) s2(80) method(asymp) `normopt'
    matrix WNS_TB=e(estimates)
    mata: st_numscalar("WNS_error",max(abs(st_matrix("WNS_TB")-st_data(.,"wns_y")'*st_matrix("WNS_TW"))))
    test_assert "TDNN `normopt' estimation consumes stable signed weights" scalar(WNS_error)<2e-12
    quietly rdistnn_predict, x0matrix(WNS_X0) sefit
    matrix WNS_prediction=r(fit)
    mata: st_numscalar("WNS_error",max(abs(st_matrix("WNS_prediction")-st_matrix("WNS_TB"))))
    test_assert "TDNN `normopt' saved prediction retains stable signed weights" scalar(WNS_error)<2e-12
}
gen double wns_g=0
gen double wns_md=0
quietly plqr_dnn wns_y, x(wns_x) x0(0) d(wns_den) md(wns_md) g(wns_g) thetaeval(0) s(80) method(asymp)
scalar WNS_plqr=el(e(moment),1,1)
mata: st_numscalar("WNS_plqr_ref",rdistnn_weights_asymp((1::100),1,80)'*((.5:-(st_data(.,"wns_y"):<=0)):*st_data(.,"wns_den")))
test_assert "PLQR moments consume the same stable asymptotic weights" abs(scalar(WNS_plqr)-scalar(WNS_plqr_ref))<2e-12
