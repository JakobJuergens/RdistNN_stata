// Frozen R-compatible distance geometry and batched observation-level weights.
quietly do "src/_rdistnn_mata.ado"
do "tests/test_framework.do"
quietly do "tests/fixtures/r_geometry_fixtures.do"
test_section "R geometry and batch weights"
foreach rule in sd mad iqr range custom {
    mata {
        GX=st_matrix("RG_X"); GX0=st_matrix("RG_X0")
        GS=rdistnn_fit_geometry(GX,"l1_scaled","`rule'",("`rule'"=="custom" ? (-2,0,4) : J(1,0,.)),1e-8)
        st_numscalar("g_scale",mreldif(GS,st_matrix("RG_scale_`rule'")))
        GD=rdistnn_distance_matrix(GX,GX0,"l1_scaled",GS)
        st_numscalar("g_l1",mreldif(GD,st_matrix("RG_l1_`rule'")))
        st_numscalar("g_rank1",mreldif(rdistnn_ranks_from_dmat(GD),st_matrix("RG_rank_l1_`rule'")))
        GD=rdistnn_distance_matrix(GX,GX0,"l2_scaled",GS)
        st_numscalar("g_l2",mreldif(GD,st_matrix("RG_l2_`rule'")))
        st_numscalar("g_rank2",mreldif(rdistnn_ranks_from_dmat(GD),st_matrix("RG_rank_l2_`rule'")))
    }
    test_assert "R `rule' fitted scales" scalar(g_scale)<1e-12
    test_assert "R `rule' scaled L1" scalar(g_l1)<1e-12
    test_assert "R `rule' scaled L2" scalar(g_l2)<1e-12
    test_assert "R `rule' scaled L1 ranks" scalar(g_rank1)==0
    test_assert "R `rule' scaled L2 ranks" scalar(g_rank2)==0
}
clear
svmat double RG_X, names(x)
gen double y=x1+2*x2
gen double den=2+_n/10
replace x1=. in 2
weights_dnn x1 x2 x3, x0matrix(RG_X0) s(3) distance(l1_scaled) scalerule(iqr)
matrix GW=r(weights)
matrix GR=r(ranks)
matrix GD=r(distances)
matrix GI=r(row_ids)
test_assert "batch sample rows and targets" rowsof(GW)==5 & colsof(GW)==2
test_assert "original row mapping survives excluded rows" GI[2,1]==3
matrix G0=RG_X0[1,1...]
weights_dnn x1 x2 x3, x0matrix(G0) s(3) distance(l1_scaled) scalerule(iqr) generate(gw)
matrix GS=r(weights)
mata: st_numscalar("g_equal",mreldif(st_matrix("GS"),st_matrix("GW")[,1]))
test_assert "single/batch weights agree" scalar(g_equal)==0
capture weights_dnn x1 x2 x3, x0matrix(RG_X0) s(3) generate(bad)
test_assert "generate rejects multiple targets" _rc==198
capture weights_dnn x1 x2 x3, x0(0 2 5) s(3) distance(l2_scaled) scalerule(custom) scalevec(1 2)
test_assert "custom wrong dimension rejected" _rc==198
capture weights_dnn x1 x2 x3, x0(0 2 5) s(3) scaleeps(0)
test_assert "nonpositive epsilon rejected" _rc==198
// Every estimator uses the same fitted geometry and rank permutations.
foreach metric in l1 l1_scaled l2_scaled mahalanobis {
    foreach rule in sd custom {
        local scaleopts "scalerule(`rule')"
        if "`rule'"=="custom" local scaleopts "`scaleopts' scalevec(-2 0 4) scaleeps(1e-8)"
        np_dnn y, xcols(x1 x2 x3) x0matrix(RG_X0) s(3) distance(`metric') `scaleopts'
        matrix GN=e(model_ranks)
        matrix GB=e(estimates)
        matrix GV=e(estimate_vcov)
        mata {
            GX=st_matrix("e(model_X)"); GY=st_matrix("e(model_signals)")
            GX0=st_matrix("RG_X0"); GG=J(0,0,.)
        }
        if inlist("`metric'","l1_scaled","l2_scaled","mahalanobis") {
            mata: GG=st_matrix("e(model_geometry)")
        }
        // Independent explicit deletion reference: preserve original row ties
        // and fixed full-sample geometry; rebuild ranks and weights each time.
        mata {
            GL=J(rows(GX),rows(GX0),.)
            for(gi=1;gi<=rows(GX);gi++) {
                gkeep=selectindex((1::rows(GX)):!=gi)
                gdx=rdistnn_distance_matrix(GX[gkeep,],GX0,"`metric'",GG)
                gr=rdistnn_ranks_from_dmat(gdx)
                gw=rdistnn_weights_from_rmat(gr,1,3,"exact")
                GL[gi,]=GY[gkeep,]'*gw
            }
            GC=GL:-mean(GL)
            st_numscalar("g_cov",mreldif(st_matrix("GV"),(rows(GX)-1)/rows(GX)*cross(GC,GC)))
        }
        test_assert "`metric' `rule' fixed-geometry deletion covariance" scalar(g_cov)<1e-11
        cate_dnn y, xcols(x1 x2 x3) x0matrix(RG_X0) s(3) pseudooutcome distance(`metric') `scaleopts'
        mata: st_numscalar("g_route",mreldif(st_matrix("GN"),st_matrix("e(model_ranks)")))
        test_assert "`metric' `rule' CATE geometry" scalar(g_route)==0
        clate_dnn y, denominator(den) pseudooutcome xcols(x1 x2 x3) x0matrix(RG_X0) s(3) distance(`metric') `scaleopts'
        mata: st_numscalar("g_route",mreldif(st_matrix("GN"),st_matrix("e(model_ranks)")))
        test_assert "`metric' `rule' CLATE geometry" scalar(g_route)==0
        // Saved model contains all training geometry; active observations may
        // be absent and target profiles may differ from the fitted profiles.
        matrix NEW0=(.7,5,8)
        mata {
            gdx=rdistnn_distance_matrix(GX,st_matrix("NEW0"),"`metric'",GG)
            gw=rdistnn_weights_from_rmat(rdistnn_ranks_from_dmat(gdx),1,3,"exact")
            st_matrix("GPred",GY'*gw)
        }
        np_dnn y, xcols(x1 x2 x3) x0matrix(RG_X0) s(3) distance(`metric') `scaleopts'
        tempfile saved
        estimates save `saved', replace
        preserve
        clear
        estimates use `saved'
        rdistnn_predict, x0matrix(NEW0) sefit
        matrix GActual=r(fit)
        mata: st_numscalar("g_predict",mreldif(st_matrix("GPred"),st_matrix("GActual")))
        restore
        test_assert "`metric' `rule' serialized frozen prediction" scalar(g_predict)<1e-12
    }
}
// Supplied-coefficient PLQR retains geometry and the matching single target.
gen double d=1+_n
gen double md=0
gen double q=0
plqr_dnn y, x(x1 x2 x3) x0(0 2 5) d(d) md(md) g(q) thetaeval(1) s(3) distance(l2_scaled) scalerule(custom) scalevec(-2 0 4) scaleeps(1e-8)
matrix GP=e(model_geometry)
test_assert "PLQR preserves fitted custom scales" GP[1,1]==2 & GP[1,2]==1e-8 & GP[1,3]==4
