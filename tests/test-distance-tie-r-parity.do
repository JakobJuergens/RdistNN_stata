// Independent subset/deletion oracles, checked against actual R/native engines
// by tools/generate_distance_tie_fixtures.R before generating this fixture.
version 17
quietly do "src/_rdistnn_mata.ado"
do "tests/test_framework.do"
quietly do "tests/fixtures/r_distance_tie_fixtures.do"
test_section "Exact distance ties: cross-language subset and deletion references"

clear
svmat double DT_X, names(dtx)
svmat double DT_y, names(dty)
svmat double DT_den, names(dtden)
svmat double DT_D, names(dtreg)
gen double dtg = 0
gen double dtmd1 = .25
gen double dtmd2 = .25

foreach metric in l1 l2 l1_scaled l2_scaled mahalanobis {
    local geometry "distance(`metric')"
    if inlist("`metric'","l1_scaled","l2_scaled") local geometry "`geometry' scalerule(custom) scalevec(2 3)"
    foreach method in exact asymp {
        foreach k in 1 2 {
            local id "`metric'_`method'`k'"
            foreach engine in mata plugin naive {
                global rdistnn_jackknife_engine `engine'
                global rdistnn_distance_engine mata
                if "`engine'"=="plugin" global rdistnn_distance_engine plugin
                quietly weights_dnn dtx1 dtx2, x0matrix(DT_X0) k(`k') s(5) method(`method') `geometry'
                local distance_route "$rdistnn_last_distance_engine"
                matrix dt_weights = r(weights)
                matrix dt_ranks = r(ranks)
                matrix dt_distances = r(distances)
                mata {
                    dt_w=st_matrix("dt_weights"); dt_d=st_matrix("dt_distances")
                    dt_refd=st_matrix("DT_`metric'_dx"); dt_same=1
                    for(dt_j=1;dt_j<=cols(dt_d);dt_j++) {
                        for(dt_i=1;dt_i<=rows(dt_d);dt_i++) {
                            dt_same=dt_same & all((dt_d[,dt_j]:==dt_d[dt_i,dt_j]):==
                                (dt_refd[,dt_j]:==dt_refd[dt_i,dt_j]))
                        }
                    }
                    st_numscalar("dt_ties",dt_same)
                    st_numscalar("dt_wr",max(abs(dt_w-st_matrix("DT_`id'_W"))))
                    st_numscalar("dt_rr",all(st_matrix("dt_ranks"):==st_matrix("DT_`metric'_ranks")))
                    st_numscalar("dt_dr",max(abs(dt_d-dt_refd)))
                }
                test_assert "`id'/`engine': actual distance engine" "`distance_route'"=="$rdistnn_distance_engine"
                test_assert "`id'/`engine': exact ranks and tie membership" scalar(dt_rr)==1 & scalar(dt_ties)==1
                test_assert "`id'/`engine': independent averaged weights and metric" scalar(dt_wr)<2e-12 & scalar(dt_dr)<2e-12

                quietly np_dnn dty1, xcols(dtx1 dtx2) x0matrix(DT_X0) k(`k') s(5) method(`method') `geometry'
                local jk_route "$rdistnn_last_jackknife_engine"
                mata: st_numscalar("dt_b",max(abs(st_matrix("e(b)")-st_matrix("DT_`id'_fit")))); st_numscalar("dt_v",max(abs(st_matrix("e(V)")-st_matrix("DT_`id'_V"))))
                test_assert "`id'/`engine': NP actual deletion engine" "`jk_route'"=="`engine'"
                test_assert "`id'/`engine': NP mean and common-target covariance" scalar(dt_b)<2e-11 & scalar(dt_v)<2e-11

                quietly rdistnn_results, type(weights)
                mata {
                    dt_data=st_matrix("r(data)"); dt_error=0
                    for(dt_j=1;dt_j<=4;dt_j++) dt_error=max((dt_error,
                        max(abs(dt_data[(dt_j-1)*12:+(1::12),5]-st_matrix("DT_`id'_W")[,dt_j]))))
                    st_numscalar("dt_reuse",dt_error)
                }
                test_assert "`id'/`engine': retained weight reconstruction" scalar(dt_reuse)<2e-12
                tempfile saved
                quietly estimates save "`saved'", replace
                preserve
                clear
                quietly estimates use "`saved'"
                quietly rdistnn_predict, x0matrix(DT_NEW) sefit
                mata: st_numscalar("dt_pb",max(abs(st_matrix("r(fit)")-st_matrix("DT_`id'_pred")))); st_numscalar("dt_pv",max(abs(st_matrix("r(vcov)")-st_matrix("DT_`id'_predV"))))
                test_assert "`id'/`engine': cold saved new-target prediction and covariance" scalar(dt_pb)<2e-11 & scalar(dt_pv)<2e-11 & _N==0
                restore

                quietly cate_dnn dty1, xcols(dtx1 dtx2) x0matrix(DT_X0) pseudooutcome k(`k') s(5) method(`method') `geometry'
                mata: st_numscalar("dt_b",max(abs(st_matrix("e(b)")-st_matrix("DT_`id'_fit")))); st_numscalar("dt_v",max(abs(st_matrix("e(V)")-st_matrix("DT_`id'_V"))))
                test_assert "`id'/`engine': CATE signal mean and covariance" scalar(dt_b)<2e-11 & scalar(dt_v)<2e-11 & "$rdistnn_last_jackknife_engine"=="`engine'"

                quietly clate_dnn dty1, xcols(dtx1 dtx2) x0matrix(DT_X0) pseudooutcome denominator(dtden1) k(`k') s(5) method(`method') `geometry'
                mata: st_numscalar("dt_b",max(abs(st_matrix("e(estimates)")-st_matrix("DT_`id'_ratio")))); st_numscalar("dt_v",max(abs(st_matrix("e(estimate_vcov)")-st_matrix("DT_`id'_ratioV"))))
                test_assert "`id'/`engine': CLATE paired components and ratio covariance" scalar(dt_b)<2e-11 & scalar(dt_v)<2e-11 & "$rdistnn_last_jackknife_engine"=="`engine'"

                quietly plqr_dnn dty1, x(dtx1 dtx2) x0(0 0) d(dtreg1 dtreg2) g(dtg) md(dtmd1 dtmd2) thetaeval(.25 -.5) k(`k') s(5) method(`method') `geometry'
                quietly rdistnn_predict, x0matrix(DT_X0) sefit
                mata: st_numscalar("dt_b",max(abs(st_matrix("r(fit)")-st_matrix("DT_`id'_moment")))); st_numscalar("dt_v",max(abs(st_matrix("r(vcov)")-st_matrix("DT_`id'_momentV"))))
                test_assert "`id'/`engine': PLQR moment vector and interleaved covariance" scalar(dt_b)<2e-11 & scalar(dt_v)<2e-11 & "$rdistnn_last_jackknife_engine"=="`engine'"
            }
        }
        local id "`metric'_tdnn_`method'"
        foreach engine in mata plugin naive {
            global rdistnn_jackknife_engine `engine'
            global rdistnn_distance_engine mata
            if "`engine'"=="plugin" global rdistnn_distance_engine plugin
            quietly np_tdnn dty1, xcols(dtx1 dtx2) x0matrix(DT_X0) s1(3) s2(6) dim(2) method(`method') `geometry'
            mata: st_numscalar("dt_b",max(abs(st_matrix("e(b)")-st_matrix("DT_`id'_fit")))); st_numscalar("dt_v",max(abs(st_matrix("e(V)")-st_matrix("DT_`id'_V"))))
            test_assert "`id'/`engine': signed TDNN mean and covariance" scalar(dt_b)<2e-11 & scalar(dt_v)<2e-11 & "$rdistnn_last_jackknife_engine"=="`engine'"
            quietly rdistnn_results, type(weights)
            mata {
                dt_data=st_matrix("r(data)"); dt_error=0
                for(dt_j=1;dt_j<=4;dt_j++) dt_error=max((dt_error,
                    max(abs(dt_data[(dt_j-1)*12:+(1::12),5]-st_matrix("DT_`id'_W")[,dt_j]))))
                st_numscalar("dt_reuse",dt_error)
            }
            test_assert "`id'/`engine': signed retained target weights" scalar(dt_reuse)<2e-12
        }
    }
}
capture macro drop rdistnn_distance_engine
capture macro drop rdistnn_jackknife_engine
