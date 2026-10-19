version 17
do "tests/test_framework.do"
quietly do "tests/fixtures/r_covariance_geometry_fixtures.do"
test_section "Actual R native/fallback geometry and covariance fixtures"
foreach case in near wide singular huge tiny offset null {
    foreach engine in mata plugin {
        global rdistnn_distance_engine `engine'
        mata {
            cgX=st_matrix("CG_`case'_X"); cgA=st_matrix("CG_`case'_A")
            cgG=rdistnn_fit_geometry(cgX,"mahalanobis","sd",J(1,0,.),1e-12)
            cgD=rdistnn_distance_matrix(cgX,cgA,"mahalanobis",cgG)
            st_numscalar("cg_distance",mreldif(cgD,st_matrix("CG_`case'_D")))
            st_numscalar("cg_ranks",all(rdistnn_ranks_from_dmat(cgD):==st_matrix("CG_`case'_ranks")))
            st_numscalar("cg_rank",cgG[1,2])
        }
        test_assert "`case'/`engine' actual projected route" "$rdistnn_last_distance_engine"=="`engine'"
        test_assert "`case'/`engine' R rank decisions" scalar(cg_ranks)==1 & scalar(cg_rank)==CG_`case'_rank[1,1]
        // Conditioning amplifies factorization rounding in the near case.
        local tol=2e-13
        if "`case'"=="near" local tol=2e-8
        test_assert "`case'/`engine' R distances" scalar(cg_distance)<`tol'
    }
}
macro drop rdistnn_distance_engine
foreach engine in mata plugin naive {
    global rdistnn_jackknife_engine `engine'
    // The largest unit-scale TDNN variance is 9.34: 2^510 would exceed
    // Stata's finite range. Keep every final variance representable here;
    // test-covariance-range separately exercises raw-sum overflow recovery.
    foreach power in 509 -510 {
        clear
        svmat double CG_X, names(cgx)
        svmat double CG_Y, names(cgy)
        gen double y=cgy1*2^`power'
        gen double den=2
        foreach method in exact asymp {
            foreach kind in np tdnn cate ratio {
                local opts "xcols(cgx1 cgx2) x0matrix(CG_A) method(`method')"
                if "`kind'"=="np" quietly np_dnn y, `opts' s(3)
                if "`kind'"=="tdnn" quietly np_tdnn y, `opts' s1(2) s2(3)
                if "`kind'"=="cate" quietly cate_dnn y, `opts' s(3) pseudooutcome
                if "`kind'"=="ratio" quietly clate_dnn y, `opts' s(3) denominator(den) pseudooutcome
                local id "`kind'_`method'"
                quietly rdistnn_results
                mata: st_numscalar("cg_report",all(st_matrix("r(fit)"):==st_matrix("e(estimates)")) & all(st_matrix("r(se_fit)"):==st_matrix("e(se_targets)")))
                test_assert "`engine'/2^`power'/`id' reporting retains numeric estimates and SEs" scalar(cg_report)==1
                mata: st_numscalar("cg_se",mreldif(sqrt(diagonal(st_matrix("e(estimate_vcov)")))'/2^`power',st_matrix("CG_`id'_se")))
                test_assert "`engine'/2^`power'/`id' R marginal covariance" scalar(cg_se)<4e-13
                if inlist("`kind'","np","tdnn") {
                    mata: st_numscalar("cg_V",mreldif(st_matrix("e(estimate_vcov)")/2^(2*`power'),st_matrix("CG_`id'_V")))
                    test_assert "`engine'/2^`power'/`id' R joint covariance" scalar(cg_V)<4e-13
                    tempfile saved
                    quietly estimates save `saved', replace
                    preserve
                    clear
                    estimates use `saved'
                    quietly rdistnn_predict, x0matrix(CG_B) sefit
                    mata: st_numscalar("cg_pred",mreldif(st_matrix("r(se_fit)")/2^`power',st_matrix("CG_`id'_pred_se")))
                    restore
                    test_assert "`engine'/2^`power'/`id' restored R covariance" scalar(cg_pred)<4e-13
                }
            }
        }
    }
}
macro drop rdistnn_jackknife_engine
