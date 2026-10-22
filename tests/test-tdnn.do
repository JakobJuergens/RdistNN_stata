// TDNN exact/asymptotic parity, signed support, centered deletion covariance.
do "tests/test_framework.do"
test_section "TDNN exact components retain canonical R signed arithmetic"
clear
set obs 7
gen double td_canonical_x=_n
gen double td_canonical_y=cond(_n==2,2^54,0)
foreach engine in mata plugin {
    global rdistnn_distance_engine `engine'
    quietly weights_tdnn td_canonical_x, x0(0) s1(3) s2(6) dim(2)
    matrix TD_canonical=r(rank_weights)
    // Portable exact component tables and separate double products give
    // -2^-54 at rank2 in both R engines; no epsilon discards this entry.
    test_assert "`engine': represented canonical rank2 sign and value" TD_canonical[2,1]==-2^-54
    quietly np_tdnn td_canonical_y, xcols(td_canonical_x) x0(0) s1(3) s2(6) dim(2) nose
    matrix TD_canonical_fit=e(estimates)
    test_assert "`engine': NP consumes canonical signed table" TD_canonical_fit[1,1]==-1
}
macro drop rdistnn_distance_engine
quietly do "tests/fixtures/r_tdnn_fixtures.do"
clear
input double x1 double x2
-1 0
0 -1
0 0
1 0
0 1
1 1
-1 -1
2 -.5
-.5 1.5
.75 -1.25
1.5 .5
-1.5 .75
end
gen double y=.7*x1-.4*x2+_n/20
matrix TD_X0=(.2,-.1 \ -.5,.75 \ 1.2,.4)
test_section "TDNN: R tables, batch estimates and full covariance"
foreach method in exact asymp {
    quietly np_tdnn y, xcols(x1 x2) x0matrix(TD_X0) s1(3) s2(6) method(`method')
    matrix TD_B=e(estimates)
    matrix TD_V=e(estimate_vcov)
    matrix TD_W=e(model_ranks)
    mata: st_numscalar("TD_error",max(abs(st_matrix("TD_B")-st_matrix("R_tdnn_`method'_fit"))))
    test_assert "`method' batch estimates match R" scalar(TD_error)<1e-11
    mata: st_numscalar("TD_error",max(abs(st_matrix("TD_V")-st_matrix("R_tdnn_`method'_vcov"))))
    test_assert "`method' full covariance matches independent R deletions" scalar(TD_error)<1e-11
    quietly np_tdnn, level(90)
    matrix TD_summary=r(summary)
    test_assert "replay exposes canonical numeric summary" rowsof(TD_summary)==3 & colsof(TD_summary)==6
    quietly weights_tdnn x1 x2, x0matrix(TD_X0) s1(3) s2(6) method(`method')
    matrix TD_weights=r(weights)
    matrix TD_table=r(rank_weights)
    mata: st_numscalar("TD_error",max(abs(st_matrix("TD_weights")-st_matrix("R_tdnn_`method'_weights"))))
    test_assert "`method' observation weights match R" scalar(TD_error)<1e-11
    mata: st_numscalar("TD_error",max(abs(st_matrix("TD_table")-st_matrix("R_tdnn_`method'_table"))))
    test_assert "`method' rank tables match R" scalar(TD_error)<1e-11
    mata: st_numscalar("TD_signed",min(st_matrix("TD_weights"))<0)
    test_assert "`method' weights retain negative support" scalar(TD_signed)==1
    quietly np_tdnn y, xcols(x1 x2) x0(.2 -.1) s1(3) s2(6) method(`method')
    test_assert "`method' scalar equals batch" abs(el(e(b),1,1)-el(TD_B,1,1))<1e-12
    global rdistnn_jackknife_engine naive
    quietly np_tdnn y, xcols(x1 x2) x0matrix(TD_X0) s1(3) s2(6) method(`method')
    mata: st_numscalar("TD_error",max(abs(st_matrix("e(estimate_vcov)")-st_matrix("TD_V"))))
    test_assert "`method' explicit deletion covariance equals prefix" scalar(TD_error)<1e-12
    global rdistnn_jackknife_engine mata
}
global rdistnn_jackknife_engine
test_section "TDNN: arbitrary signed table covariance and support centering"
quietly run "src/_rdistnn_tdnn_mata.ado"
mata {
    TD_R=((1::8),(8::1)); TD_y=(1e12:+(0::6)/1024)\1e300
    TD_wn=(.7\.4\-.1\0\0\0\0\0)
    TD_wd=(.7\.4\-.1\0\0\0\0)
    TD_fast=rdistnn_loo_from_table(TD_R[,1],TD_y,TD_wd,1,"mata")
    TD_slow=rdistnn_loo_from_table(TD_R[,1],TD_y,TD_wd,1,"naive")
    st_numscalar("TD_error",max(abs(TD_fast-TD_slow)))
    st_numscalar("TD_finite",!hasmissing(TD_fast))
    TD_const=rdistnn_loo_from_table(TD_R,J(8,1,1e200),TD_wd,1,"mata")
    st_numscalar("TD_constant",max(abs(TD_const)))
}
test_assert "arbitrary signed table prefix equals independent deletions" scalar(TD_error)<1e-12
test_assert "unused huge outcome does not overflow centered deletions" scalar(TD_finite)==1
test_assert "large constant outcomes yield exact zero covariance inputs" scalar(TD_constant)==0
test_section "TDNN: common signed-table engine dispatch"
foreach engine in mata naive plugin {
    if "`engine'"=="plugin" & "$rdistnn_plugin_ok"!="1" continue
    global rdistnn_jackknife_engine `engine'
    foreach method in exact asymp {
        quietly np_tdnn y, xcols(x1 x2) x0matrix(TD_X0) s1(3) s2(6) method(`method')
        matrix TD_engineV=e(estimate_vcov)
        mata: st_numscalar("TD_error",max(abs(st_matrix("TD_engineV")-st_matrix("R_tdnn_`method'_vcov"))))
        test_assert "`engine' `method' signed-table covariance matches R" scalar(TD_error)<1e-11
    }
}
capture macro drop rdistnn_jackknife_engine
gen double td_constant=1e200
quietly np_tdnn td_constant, xcols(x1 x2) x0matrix(TD_X0) s1(3) s2(6)
mata: st_numscalar("TD_constant",max(abs(st_matrix("e(estimate_vcov)"))))
test_assert "constant-outcome TDNN full covariance is exactly zero" scalar(TD_constant)==0
quietly np_tdnn, level(90)
matrix TD_summary=r(summary)
test_assert "zero standard error gives collapsed confidence bounds" el(TD_summary,1,5)==el(TD_summary,1,1) & el(TD_summary,1,6)==el(TD_summary,1,1)
quietly weights_tdnn x1 x2, x0matrix(TD_X0) s1(3) s2(6) dim(4) method(asymp) nonormalize
local TD_normalized=r(normalize)
matrix TD_custom=r(rank_weights)
test_assert "asymptotic nonormalize is retained" `TD_normalized'==0
mata: st_numscalar("TD_error",max(abs(st_matrix("TD_custom")-st_matrix("R_tdnn_asymp_dim4_nonorm"))))
test_assert "explicit bias dimension and nonormalize match R" scalar(TD_error)<1e-12
test_section "TDNN: boundaries, ties and signed diagnostics"
quietly np_tdnn y, xcols(x1 x2) x0(0 0) s1(3) s2(12) nose
test_assert "nose permits s2=n" e(s2)==12 & "`e(se_method)'"=="none"
test_assert "nose preserves missing uncertainty" missing(el(e(se_targets),1,1))
capture np_tdnn y, xcols(x1 x2) x0(0 0) s1(3) s2(12)
test_assert "jackknife rejects s2=n" _rc==198
capture np_tdnn y, xcols(x1 x2) x0(0 0) s1(6) s2(3)
test_assert "reject reversed scales" _rc==198
capture weights_tdnn x1 x2, x0(0) s1(3) s2(6)
test_assert "reject wrong target dimensions" _rc==198
capture weights_tdnn x1 x2, x0(0 0) s1(3) s2(6) dim(-1)
test_assert "reject negative bias dimension" _rc==198
capture weights_tdnn x1 x2, x0(0 0) s1(3) s2(6) dim(0)
test_assert "reject zero bias dimension" _rc==198
quietly weights_tdnn x1 x2, x0(0 0) s1(3) s2(6) generate(td_w)
matrix TD_r=r(ranks)
test_assert "stable tied distance order" el(TD_r,1,1)<el(TD_r,2,1)
quietly np_tdnn y, xcols(x1 x2) x0matrix(TD_X0) s1(3) s2(6)
quietly np_tdnn_diagnostics, targets(2 1)
matrix TD_d=r(diagnostics)
test_assert "signed diagnostics retain requested order" el(TD_d,1,1)==2 & el(TD_d,2,1)==1
test_assert "signed diagnostics report negative mass" el(TD_d,1,6)>0
quietly np_tdnn_plot, type(weights) nodraw replace
matrix TD_plot=r(plot_data)
test_assert "signed weight plot exposes numeric data" colsof(TD_plot)==3
test_section "TDNN: observation mapping and signed concentration state"
preserve
replace x1=. in 4
quietly weights_tdnn x1 x2 if _n!=2 in 2/11, x0(0 0) s1(2) s2(4)
matrix TD_rows=r(row_ids)
matrix TD_weights=r(weights)
mata: st_numscalar("TD_rows_ok",all(st_matrix("TD_rows"):==(3\5\6\7\8\9\10\11)))
test_assert "weight rows map to original observations after exclusions" scalar(TD_rows_ok)==1 & rowsof(TD_rows)==rowsof(TD_weights)
restore
quietly np_tdnn y, xcols(x1 x2) x0matrix(TD_X0) s1(3) s2(6)
matrix TD_concentration=e(neff_targets)
matrix TD_diag=e(diagnostics)
mata: st_numscalar("TD_error",max(abs(st_matrix("TD_concentration")'-st_matrix("TD_diag")[,5])))
test_assert "signed concentration state agrees with diagnostics" scalar(TD_error)==0 & !missing(e(neff_min),e(neff_max))
quietly np_tdnn y, xcols(x1 x2) x0(0 0) s1(3) s2(6)
test_assert "scalar signed concentration is retained for replay" !missing(e(neff)) & e(neff)==el(e(diagnostics),1,5)
