// Actual R native/fallback references for weights, extreme geometry and reuse.
version 17
quietly do "src/_rdistnn_mata.ado"
quietly do "src/_rdistnn_tdnn_mata.ado"
do "tests/test_framework.do"
quietly do "tests/fixtures/r_stability_fixtures.do"
test_section "Numerical stability: R tail tables and close-order coefficients"
foreach k in 1 3 80 {
    mata: st_numscalar("st_error",max(abs(rdistnn_weights_asymp((1::100),`k',80):/st_matrix("ST_W`k'"):-1)))
    test_assert "k`k': represented tail weights agree relatively with R" scalar(st_error)<2e-10
}
mata {
    st_c=rdistnn_tdnn_coefficients(1000,1001,50)\rdistnn_tdnn_coefficients(10000,10001,100)
    st_numscalar("st_error",max(abs(st_c:/st_matrix("ST_COEF"):-1)))
}
test_assert "close-order coefficients agree with R" scalar(st_error)<2e-14
foreach method in exact asymp {
    // Coefficients near 4000 amplify existing lncomb/lchoose rounding and can
    // cross the retained normalization threshold. Bound continuous error on
    // the table's scale and check an independent 90-digit exact-rule oracle.
    mata {
        st_w=rdistnn_tdnn_table(100,79,80,100,"`method'",1)
        st_ref=st_matrix("ST_TD_`method'")
        st_numscalar("st_error",max(abs(st_w-st_ref)))
        st_numscalar("st_bound",("`method'"=="exact") ?
            1e-9*max((1,max(abs(st_ref)))) : 2e-10)
        st_numscalar("st_support",all((st_w:!=0):==(st_ref:!=0)))
        st_numscalar("st_sign",all(sign(st_w):==sign(st_ref)))
    }
    test_assert "`method': normalized signed table agrees with R" scalar(st_error)<scalar(st_bound)
    test_assert "`method': signed support agrees exactly with R" scalar(st_support)==1 & scalar(st_sign)==1
    if "`method'"=="exact" {
        mata: st_numscalar("st_oracle",max(abs(st_w-st_matrix("ST_TD_exact_oracle")))/max(abs(st_matrix("ST_TD_exact_oracle"))))
        test_assert "exact: independent high-precision signed table" scalar(st_oracle)<1e-9
    }
}
test_section "Numerical stability: fitted R parity and empty-data saved prediction"
global rdistnn_jackknife_engine mata
foreach size in huge tiny {
    clear
    svmat double ST_`size'_X, names(stx)
    svmat double ST_`size'_Y, names(sty)
    local h "1e200"
    if "`size'"=="tiny" local h "1e-200"
    foreach metric in l2 l2_scaled {
        local geometry "distance(`metric')"
        if "`metric'"=="l2_scaled" local geometry "`geometry' scalerule(custom) scalevec(2 4)"
        foreach method in exact asymp {
            local id "`size'_`metric'_`method'"
            foreach engine in mata plugin {
                global rdistnn_distance_engine `engine'
                quietly np_dnn sty1, xcols(stx1 stx2) x0matrix(ST_`size'_A) k(2) s(4) method(`method') `geometry'
                local route "$rdistnn_last_distance_engine"
                mata {
                    st_numscalar("st_fit",max(abs(st_matrix("e(b)")-st_matrix("ST_`id'_fit"))))
                    st_numscalar("st_V",max(abs(st_matrix("e(V)")-st_matrix("ST_`id'_V"))))
                    st_numscalar("st_D",max(abs((st_matrix("e(model_distances)")-st_matrix("ST_`id'_D")):/`h')))
                    st_numscalar("st_ranks",all(st_matrix("e(model_ranks)"):==st_matrix("ST_`id'_ranks")))
                }
                test_assert "`id'/`engine': real route and exact R ranks" "`route'"=="`engine'" & scalar(st_ranks)==1
                test_assert "`id'/`engine': R distances, estimate and joint covariance" scalar(st_D)<2e-13 & scalar(st_fit)<2e-11 & scalar(st_V)<2e-11
                tempfile saved
                quietly estimates save "`saved'", replace
                preserve
                clear
                quietly estimates use "`saved'"
                quietly rdistnn_predict, x0matrix(ST_`size'_B) sefit
                local route "$rdistnn_last_distance_engine"
                mata {
                    st_numscalar("st_p",max(abs(st_matrix("r(fit)")-st_matrix("ST_`id'_pred"))))
                    st_numscalar("st_s",max(abs(st_matrix("r(se_fit)")-st_matrix("ST_`id'_pse"))))
                }
                test_assert "`id'/`engine': restored R prediction and uncertainty" scalar(st_p)<2e-11 & scalar(st_s)<2e-11 & _N==0 & "`route'"=="`engine'"
                restore
            }
        }
    }
}
capture macro drop rdistnn_distance_engine
capture macro drop rdistnn_last_distance_engine
capture macro drop rdistnn_jackknife_engine
