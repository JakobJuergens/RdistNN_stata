// Independent analytic scales and R-backed neighbor decisions after fitting.
version 17
quietly do "src/_rdistnn_mata.ado"
do "tests/test_framework.do"
quietly do "tests/fixtures/r_geometry_fixtures.do"
test_section "Fitted geometry: guarded R quantile arithmetic"
mata {
    gs_QX=st_matrix("RG_IQR_X"); gs_QT=st_matrix("RG_IQR_X0")
    gs_G=rdistnn_fit_geometry(gs_QX,"l2_scaled","iqr",J(1,0,.),1e-12)
    st_numscalar("gs_quantiles",rdistnn_type7((-8e307\8e307),0)==-8e307 &
        rdistnn_type7((-8e307\8e307),.5)==0 &
        rdistnn_type7((-8e307\8e307),1)==8e307 &
        rdistnn_type7((2\2),.25)==2 & rdistnn_type7(3,.75)==3)
}
test_assert "type7 endpoints, midpoint, repeats and singleton remain finite" scalar(gs_quantiles)==1
foreach engine in mata plugin {
    global rdistnn_distance_engine `engine'
    mata {
        gs_D=rdistnn_distance_matrix(gs_QX,gs_QT,"l2_scaled",gs_G)
        st_numscalar("gs_ranks",all(rdistnn_ranks_from_dmat(gs_D):==st_matrix("RG_IQR_rank")))
        st_numscalar("gs_distances",mreldif(gs_D,st_matrix("RG_IQR_D")))
    }
    test_assert "`engine': IQR near-tie ranks equal R exactly" scalar(gs_ranks)==1
    test_assert "`engine': IQR distances equal R" scalar(gs_distances)<2e-15
    test_assert "`engine': IQR actual distance route" "$rdistnn_last_distance_engine"=="`engine'"
    clear
    svmat double RG_IQR_X, names(gsx)
    gen double gsy=cond(_n==1,100,0)
    quietly np_dnn gsy, xcols(gsx1 gsx2) x0matrix(RG_IQR_X0) s(10) nose distance(l2_scaled) scalerule(iqr)
    test_assert "`engine': IQR near-tie actual estimate equals R" el(e(b),1,1)==100
}

test_section "Fitted geometry: SD recovery keeps R near-tie neighbors"
foreach engine in mata plugin {
    global rdistnn_distance_engine `engine'
    mata {
        gs_X=st_matrix("RG_SD_X"); gs_T=st_matrix("RG_SD_X0")
        gs_G=rdistnn_fit_geometry(gs_X,"l2_scaled","sd",J(1,0,.),1e-12)
        gs_D=rdistnn_distance_matrix(gs_X,gs_T,"l2_scaled",gs_G)
        st_numscalar("gs_sd_rank",all(rdistnn_ranks_from_dmat(gs_D):==st_matrix("RG_SD_rank")))
        st_numscalar("gs_sd_scale",all((1:/gs_G):==st_matrix("RG_SD_inv")))
    }
    test_assert "`engine': recovered SD inverse scales equal R exactly" scalar(gs_sd_scale)==1
    test_assert "`engine': recovered SD ranks equal R exactly" scalar(gs_sd_rank)==1
    clear
    svmat double RG_SD_X, names(gsx)
    gen double gsy=cond(_n==2,100,0)
    quietly np_dnn gsy, xcols(gsx1 gsx2) x0matrix(RG_SD_X0) s(8) nose distance(l2_scaled) scalerule(sd)
    test_assert "`engine': recovered SD near-tie estimate equals R" el(e(b),1,1)==100
}

test_section "Fitted geometry: finite SD after intermediate overflow or underflow"
foreach h in 1e154 1e200 1e-200 {
    local gs_eps=cond("`h'"=="1e-200",1e-250,1e-12)
    mata {
        gs_h=strtoreal("`h'")
        gs_X=(-1,0\-1,1\1,0\1,1):*(gs_h,1)
        gs_G=rdistnn_fit_geometry(gs_X,"l2_scaled","sd",J(1,0,.),1e-250)
        st_numscalar("gs_sd",abs(gs_G[1]/gs_h-sqrt(4/3)))
        st_matrix("GS_X",gs_X)
    }
    test_assert "`h': analytic sample SD survives intermediate range" scalar(gs_sd)<5e-15
    foreach engine in mata plugin {
        global rdistnn_distance_engine `engine'
        clear
        svmat double GS_X, names(gsx)
        gen double gsy=cond(_n==2,100,0)
        quietly np_dnn gsy, xcols(gsx1 gsx2) x0(0 1) s(4) nose distance(l2_scaled) scalerule(sd) scaleeps(`gs_eps')
        test_assert "`engine'/`h': SD actual estimate equals analytic/R reference" el(e(b),1,1)==100
        tempfile gs_saved
        quietly estimates save "`gs_saved'", replace
        preserve
        clear
        quietly estimates use "`gs_saved'"
        matrix GS_NEW=(0,1\0,1)
        quietly rdistnn_predict, x0matrix(GS_NEW)
        matrix GS_P=r(fit)
        test_assert "`engine'/`h': saved fitted SD survives duplicate targets" GS_P[1,1]==100 & GS_P[1,2]==100
        restore
    }
}
capture mata: gs_G=rdistnn_fit_geometry((-6e307\6e307),"l2_scaled","range",J(1,0,.),1e-12)
local gs_rc=_rc
test_assert "unrepresentable fitted range errors before epsilon floor" `gs_rc'==498
mata {
    gs_G=rdistnn_fit_geometry((-6e307\6e307),"l2_scaled","mad",J(1,0,.),1e-12)
    st_numscalar("gs_mad",abs(gs_G[1]/6e307-1.4826))
    gs_G=rdistnn_fit_geometry((2\2),"l2_scaled","sd",J(1,0,.),1e-12)
    st_numscalar("gs_constant",gs_G[1]==1e-12)
}
test_assert "MAD interpolation retains finite analytic scale" scalar(gs_mad)<5e-15
test_assert "finite constant scale still uses epsilon" scalar(gs_constant)==1
mata {
    // Decimal literals keep the binary powers identical to the R inputs.
    gs_unit=4.6719391924451279e+195
    gs_offset=5.2601359015483735e+210:+(-1\-1\1\1):*gs_unit
    gs_G=rdistnn_fit_geometry(gs_offset,"l2_scaled","sd",J(1,0,.),1e-12)
    st_numscalar("gs_offset_sd",abs(gs_G[1]/gs_unit-sqrt(4/3)))
}
test_assert "SD range recovery preserves small variation around a huge offset" scalar(gs_offset_sd)<5e-15

test_section "Scaled L1: shared boundary, actual engines and exact neighbor order"
foreach engine in mata plugin {
    global rdistnn_distance_engine `engine'
    capture mata: gs_L=rdistnn_distance_matrix((8.000000000000001e307\8e307),(-7.9e307\-7.9e307),"l1_scaled",10)
    local gs_rc=_rc
    test_assert "`engine': finite scaled L1 accepts exceptional raw difference" `gs_rc'==0
    if `gs_rc'==0 {
        mata {
            st_numscalar("gs_l1",all(gs_L:==(1.5900000000000001e307\1.59e307)))
            st_numscalar("gs_l1_ranks",all(rdistnn_ranks_from_dmat(gs_L):==(2\1)))
            gs_single=rdistnn_distance_matrix((8.000000000000001e307\8e307),-7.9e307,"l1_scaled",10)
            st_numscalar("gs_l1_single",all(gs_single:==gs_L[,1]))
        }
        test_assert "`engine': scaled L1 exact R distances" scalar(gs_l1)==1
        test_assert "`engine': scaled L1 exact R ranks" scalar(gs_l1_ranks)==1
        test_assert "`engine': scaled L1 scalar equals batch" scalar(gs_l1_single)==1
        test_assert "`engine': scaled L1 actual route" "$rdistnn_last_distance_engine"=="`engine'"
    }
    capture mata: gs_L=rdistnn_distance_matrix((8e307,8e307),(0,0),"l1",J(0,0,.))
    local gs_rc=_rc
    test_assert "`engine': unrepresentable final L1 still errors" `gs_rc'==498
}
macro drop rdistnn_distance_engine
