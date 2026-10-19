version 17
do "tests/test_framework.do"
quietly do "tests/fixtures/r_accumulation_fixtures.do"
test_section "actual R native/interpreted accumulation fixtures"
forvalues i=1/32 {
    mata: st_numscalar("ac_sum",rdistnn_sum(st_matrix("AC_sum`i'")))
    test_assert "binary signed reduction `i' matches independent unit residual" scalar(ac_sum)==1
}
foreach method in exact asymp {
    foreach kind in dnn tdnn {
        local id "`kind'_`method'"
        foreach engine in mata plugin naive {
            mata {
                acR=st_matrix("AC_ranks");acZ=st_matrix("AC_signal")
                acTable=st_matrix("AC_`id'_table");acFull=st_matrix("AC_`id'_full")
                acL=rdistnn_loo_from_table(acR,acZ,acTable,0,"`engine'")
                acW=J(rows(acR),cols(acR),0)
                for(acj=1;acj<=cols(acR);acj++)acW[,acj]=acFull[acR[,acj]]
                st_numscalar("ac_point",max(abs(rdistnn_weighted_signals(acZ,acW)-st_matrix("AC_`id'_theta"))))
                st_numscalar("ac_loo",max(abs(acL-st_matrix("AC_`id'_loo"))))
            }
            test_assert "`id'/`engine': actual R paired weighted means" scalar(ac_point)==0
            // Explicit deletion and grouped prefix reductions can round in a
            // different order. No rank/support/status tolerance is applied.
            test_assert "`id'/`engine': actual R paired deletion values" scalar(ac_loo)<2e-15
        }
    }
}

test_section "R estimator means, covariance and restored prediction"
foreach engine in mata plugin naive {
    global rdistnn_jackknife_engine `engine'
    foreach method in exact asymp {
        foreach kind in np cate ratio tdnn {
            clear
            svmat double AC_X, names(acx)
            svmat double AC_Y, names(acy)
            svmat double AC_den, names(acden)
            local opts "xcols(acx1 acx2) x0matrix(AC_targets) method(`method')"
            if "`kind'"=="np" quietly np_dnn acy1, `opts' s(7)
            if "`kind'"=="cate" quietly cate_dnn acy1, `opts' s(7) pseudooutcome
            if "`kind'"=="ratio" quietly clate_dnn acy1, `opts' s(7) denominator(acden1) pseudooutcome
            if "`kind'"=="tdnn" quietly np_tdnn acy1, `opts' s1(4) s2(7)
            local id "`kind'_`method'"
            mata: st_numscalar("ac_point",mreldif(st_matrix("e(estimates)"),st_matrix("AC_`id'_point")))
            mata: st_numscalar("ac_se",mreldif(st_matrix("e(se_targets)"),st_matrix("AC_`id'_se")))
            test_assert "`engine'/`id': actual R point estimates" scalar(ac_point)<3e-13
            test_assert "`engine'/`id': actual R standard errors" scalar(ac_se)<3e-13
            if inlist("`kind'","np","tdnn") {
                mata: st_numscalar("ac_cov",mreldif(st_matrix("e(estimate_vcov)"),st_matrix("AC_`id'_vcov")))
                test_assert "`engine'/`id': actual R joint covariance" scalar(ac_cov)<3e-13
                tempfile saved
                quietly estimates save `saved', replace
                preserve
                clear
                estimates use `saved'
                quietly rdistnn_predict, x0matrix(AC_predict) sefit
                mata: st_numscalar("ac_pred",mreldif(st_matrix("r(fit)"),st_matrix("AC_`id'_pred")))
                mata: st_numscalar("ac_predse",mreldif(st_matrix("r(se_fit)"),st_matrix("AC_`id'_pred_se")))
                restore
                test_assert "`engine'/`id': restored R prediction with duplicate targets" scalar(ac_pred)<3e-13
                test_assert "`engine'/`id': restored R prediction uncertainty" scalar(ac_predse)<3e-13
            }
        }
        clear
        svmat double AC_X, names(acx)
        svmat double AC_Y, names(acy)
        svmat double AC_D, names(acd)
        gen double acg=0
        gen double acmd=0
        mata: st_matrix("AC_one_target",st_matrix("AC_targets")[1,])
        local point1=AC_one_target[1,1]
        local point2=AC_one_target[1,2]
        quietly plqr_dnn acy1, x(acx1 acx2) x0(`point1' `point2') ///
            d(acd1) md(acmd) g(acg) thetaeval(0) s(7) method(`method')
        mata: st_numscalar("ac_moment",mreldif(st_matrix("e(moment)"),st_matrix("AC_plqr_`method'_point")))
        mata: st_numscalar("ac_momentse",mreldif(st_matrix("e(se_targets)"),st_matrix("AC_plqr_`method'_se")))
        test_assert "`engine'/`method': supplied-coefficient R PLQR moments" scalar(ac_moment)<3e-13
        test_assert "`engine'/`method': supplied-coefficient R PLQR moment SEs" scalar(ac_momentse)<3e-13
    }
}
macro drop rdistnn_jackknife_engine
