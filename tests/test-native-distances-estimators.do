// Native vs Mata end-to-end geometry for every estimator and numeric snapshots.
version 17
quietly do "src/_rdistnn_mata.ado"
do "tests/test_framework.do"
test_section "native distances: all estimators, tuning, joint covariance and saved geometry"
matrix ND_API=(73)
capture plugin call rdistnn_plugin, distance_api ND_API
local native_rc=_rc
test_assert "estimator parity requires real native distance capability" `native_rc'==0 & ND_API[1,1]==1
if `native_rc'!=0 exit 0
// Keep jackknife engine fixed: differences below isolate the distance route.
global rdistnn_jackknife_engine mata
global rdistnn_use_plugin 1
clear
set obs 24
gen double nd_x1=(_n-12)/5
gen double nd_x2=sin(_n)
gen double nd_x3=5
gen double nd_y=sin(nd_x1)+nd_x2/3+_n/50
gen double nd_den=1+cos(nd_x1)/5
gen double nd_d1=.3*nd_x1+sin(_n/3)
gen double nd_d2=-.2*nd_x1+cos(_n/3)
gen double nd_z=0
// Complete sample and original row positions must survive both exclusions.
replace nd_x1=. in 2
matrix ND_A=(-.5,0,5 \ 0,.2,5 \ -.5,0,5)
matrix ND_B=(.7,-.1,5 \ -.2,.1,5 \ .7,-.1,5)
foreach metric in l1 l2 l1_scaled l2_scaled {
    local rules "sd"
    if strpos("`metric'","scaled") local rules "sd mad iqr range custom"
    foreach rule of local rules {
        local geo "distance(`metric') scalerule(`rule') scaleeps(1e-8)"
        if "`rule'"=="custom" local geo "`geo' scalevec(-2 0 4)"
        foreach method in exact asymp {
            foreach family in np cate clate tdnn plqr {
                local cmd "`family'_dnn"
                local base "xcols(nd_x1 nd_x2 nd_x3) x0matrix(ND_A) k(2) `geo' method(`method')"
                local grid "tuningmode(grid) sgrid(3 6 4 5) sdefault(4)"
                local call "`cmd' nd_y if _n!=5, `base' `grid'"
                if "`family'"=="cate" local call "`call' pseudooutcome"
                if "`family'"=="clate" local call "`call' pseudooutcome denominator(nd_den)"
                if "`family'"=="tdnn" local call "np_tdnn nd_y if _n!=5, xcols(nd_x1 nd_x2 nd_x3) x0matrix(ND_A) s1(3) s2(6) `geo' method(`method')"
                if "`family'"=="plqr" local call "plqr_dnn nd_y if _n!=5, x(nd_x1 nd_x2 nd_x3) x0(-.5 0 5) d(nd_d1 nd_d2) md(nd_z nd_z) g(nd_z) thetaeval(.3 -.2) k(2) s(4) `geo' method(`method')"
                global rdistnn_distance_engine mata
                quietly `call'
                local observed "$rdistnn_last_distance_engine"
                matrix ND_BREF=e(b)
                matrix ND_VREF=e(V)
                matrix ND_DREF=e(model_distances)
                matrix ND_RREF=e(model_ranks)
                matrix ND_IREF=e(model_row_ids)
                if strpos("`metric'","scaled") matrix ND_GREF=e(model_geometry)
                if inlist("`family'","np","cate","clate") {
                    matrix ND_TREF=e(tuning_estimates)
                    matrix ND_TSREF=e(tuning_se)
                    quietly rdistnn_candidate, candidate(2)
                    matrix ND_CREF=r(fit)
                    matrix ND_CVREF=r(vcov)
                }
                test_assert "`family'/`metric'/`rule'/`method': reference actually uses Mata" "`observed'"=="mata"
                // Distinct new targets, including a duplicate, exercise frozen
                // training scales instead of recomputing from query profiles.
                quietly rdistnn_predict, x0matrix(ND_B) sefit
                matrix ND_PREF=r(fit)
                matrix ND_PVREF=r(vcov)
                global rdistnn_distance_engine plugin
                quietly `call'
                local observed "$rdistnn_last_distance_engine"
                mata {
                    st_numscalar("nd_fit",max(abs(st_matrix("e(b)")-st_matrix("ND_BREF"))))
                    st_numscalar("nd_cov",max(abs(st_matrix("e(V)")-st_matrix("ND_VREF"))))
                    st_numscalar("nd_dist",max(abs(st_matrix("e(model_distances)")-st_matrix("ND_DREF"))):/(1+max(abs(st_matrix("ND_DREF")))))
                    st_numscalar("nd_rank",max(abs(st_matrix("e(model_ranks)")-st_matrix("ND_RREF"))))
                    st_numscalar("nd_ids",max(abs(st_matrix("e(model_row_ids)")-st_matrix("ND_IREF"))))
                }
                test_assert "`family'/`metric'/`rule'/`method': fit actually uses native geometry" "`observed'"=="plugin"
                test_assert "`family'/`metric'/`rule'/`method': fit and full covariance match Mata" scalar(nd_fit)<1e-11 & scalar(nd_cov)<1e-11
                test_assert "`family'/`metric'/`rule'/`method': distances/ranks and complete-sample mapping match" scalar(nd_dist)<1e-14 & scalar(nd_rank)==0 & scalar(nd_ids)==0 & ND_IREF[2,1]==3
                if strpos("`metric'","scaled") {
                    mata: st_numscalar("nd_geometry",max(abs(st_matrix("e(model_geometry)")-st_matrix("ND_GREF"))))
                    test_assert "`family'/`metric'/`rule'/`method': fitted geometry retained" scalar(nd_geometry)==0
                }
                if inlist("`family'","np","cate","clate") {
                    mata: st_numscalar("nd_grid",max(abs(st_matrix("e(tuning_estimates)")-st_matrix("ND_TREF")))+max(abs(st_matrix("e(tuning_se)")-st_matrix("ND_TSREF"))))
                    test_assert "`family'/`metric'/`rule'/`method': candidate estimates and uncertainty match" scalar(nd_grid)<1e-11
                    quietly rdistnn_candidate, candidate(2)
                    mata: st_numscalar("nd_grid",max(abs(st_matrix("r(fit)")-st_matrix("ND_CREF")))+max(abs(st_matrix("r(vcov)")-st_matrix("ND_CVREF"))))
                    test_assert "`family'/`metric'/`rule'/`method': candidate numeric reuse matches" scalar(nd_grid)<1e-11
                }
                tempfile nd_saved
                quietly estimates save "`nd_saved'", replace
                preserve
                clear
                quietly estimates use "`nd_saved'"
                quietly rdistnn_predict, x0matrix(ND_B) sefit
                local observed "$rdistnn_last_distance_engine"
                mata {
                    st_numscalar("nd_predict",max(abs(st_matrix("r(fit)")-st_matrix("ND_PREF"))))
                    st_numscalar("nd_predictv",max(abs(st_matrix("r(vcov)")-st_matrix("ND_PVREF"))))
                }
                test_assert "`family'/`metric'/`rule'/`method': restored prediction actually uses native route" "`observed'"=="plugin"
                test_assert "`family'/`metric'/`rule'/`method': empty-data restored predictions and joint covariance match" scalar(nd_predict)<1e-11 & scalar(nd_predictv)<1e-11 & _N==0
                if strpos("`metric'","scaled") {
                    mata: st_numscalar("nd_geometry",max(abs(st_matrix("e(model_geometry)")-st_matrix("ND_GREF"))))
                    test_assert "`family'/`metric'/`rule'/`method': query/empty data never refits scales" scalar(nd_geometry)==0
                }
                restore
            }
        }
    }
}
capture macro drop rdistnn_distance_engine
capture macro drop rdistnn_last_distance_engine
capture macro drop rdistnn_jackknife_engine
capture macro drop rdistnn_use_plugin
