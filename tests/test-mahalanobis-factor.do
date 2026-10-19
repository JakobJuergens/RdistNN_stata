// The true nearest neighbor is known in latent orthogonal coordinates.
do "tests/test_framework.do"
test_section "Factored Mahalanobis geometry"
mata {
    mpower=-600
    mvalue=rdistnn_power2(mpower)
    st_numscalar("mahal_power",mpower==-600 & mvalue>0)
}
test_assert "binary power preserves caller exponent" scalar(mahal_power)==1
mata {
    MPX=(2^54,1,-2^54 \ 1,2,3)
    MPB=(1,0,0 \ 1,0,0 \ 1,0,0)
    MPG=(1,1,0 \ 0,0,0 \ MPB)
    MPE=(1,0,0 \ 6,0,0)
    st_numscalar("mahal_product",all(rdistnn_mahal_project(MPX,MPG):==MPE))
}
test_assert "whitening projection retains unit cancellation residual" scalar(mahal_product)==1
mata {
    mx=(-2\-1\0\1\2); mz=(1\-1\0\-1\1)
    MX=(mx,mx:+2^-22*mz)
    ma=(-17/8-2^-10 \ -17/8+2^-10)
    MT=(ma,ma:-2^-22/8)
    MG=rdistnn_fit_geometry(MX,"mahalanobis","sd",J(1,0,.),1e-12)
    ME=J(5,2,.)
    for(mj=1;mj<=2;mj++) ME[,mj]=sqrt((2/5)*(mx:-ma[mj]):^2+(mz:+1/8):^2)
    st_matrix("MX",MX); st_matrix("MT",MT)
    st_numscalar("mahal_rank",MG[1,2])
}
test_assert "resolvable near-collinear design keeps full rank" scalar(mahal_rank)==2
foreach engine in mata plugin {
    global rdistnn_distance_engine `engine'
    mata {
        MD=rdistnn_distance_matrix(MX,MT,"mahalanobis",MG)
        st_numscalar("mahal_value",mreldif(MD,ME))
        MR=rdistnn_ranks_from_dmat(MD)
        st_numscalar("mahal_nearest",MR[1,1]==1 & MR[2,2]==1)
        st_numscalar("mahal_subset",mreldif(rdistnn_distance_matrix(MX[(3\1),],MT[(2\1\2),],"mahalanobis",MG),MD[(3\1),(2\1\2)]))
        MPX=rdistnn_mahal_project(MX,MG)
        MPT=rdistnn_mahal_project(MT,MG)
        mshape=1
        for(mi=1;mi<=rows(MX);mi++) mshape=mshape & all(rdistnn_mahal_project(MX[mi,],MG):==MPX[mi,])
        for(mj=1;mj<=rows(MT);mj++) {
            mshape=mshape & all(rdistnn_mahal_project(MT[mj,],MG):==MPT[mj,])
            mshape=mshape & all(rdistnn_distance_matrix(MX,MT[mj,],"mahalanobis",MG):==MD[,mj])
        }
        st_numscalar("mahal_shape",mshape)
    }
    test_assert "`engine' analytic distances" scalar(mahal_value)<2e-8
    test_assert "`engine' exact nearest rows" scalar(mahal_nearest)==1
    test_assert "`engine' used transformed L2 engine" "$rdistnn_last_distance_engine"=="`engine'"
    test_assert "`engine' frozen subset and target order" scalar(mahal_subset)==0
    test_assert "`engine' scalar/batch projection and distance identity" scalar(mahal_shape)==1
    clear
    svmat double MX, names(x)
    gen double y=cond(_n==2,100,0)
    np_dnn y, xcols(x1 x2) x0matrix(MT) s(5) distance(mahalanobis) nose
    matrix MF=e(estimates)
    test_assert "`engine' NP nearest prediction" MF[1,1]==0 & MF[1,2]==100
    tempfile saved
    estimates save `saved', replace
    clear
    estimates use `saved'
    matrix MNEW=MT[2,1...] \ MT[1,1...] \ MT[2,1...]
    rdistnn_predict, x0matrix(MNEW)
    matrix MF=r(fit)
    test_assert "`engine' serialized frozen prediction" MF[1,1]==100 & MF[1,2]==0 & MF[1,3]==100
}
macro drop rdistnn_distance_engine
mata {
    MB=(2,1\2,-1\-2,1\-2,-1\0,0)
    ME=sqrt(MB[,1]:^2/4+MB[,2]:^2)
    for(ms=-600;ms<=600;ms=ms+1200) {
        MX=MB*2^ms
        MG=rdistnn_fit_geometry(MX,"mahalanobis","sd",J(1,0,.),1e-12)
        st_numscalar("mahal_extreme_"+strofreal(ms==600),mreldif(rdistnn_distance_matrix(MX,(0,0),"mahalanobis",MG),ME))
    }
    MX=(mx,mx)
    MG=rdistnn_fit_geometry(MX,"mahalanobis","sd",J(1,0,.),1e-12)
    st_numscalar("mahal_singular",MG[1,2]==1 & mreldif(rdistnn_distance_matrix(MX,(1,-1),"mahalanobis",MG),abs(mx)/sqrt(2.5))<2e-14)
    MG=rdistnn_fit_geometry(J(3,2,7),"mahalanobis","sd",J(1,0,.),1e-12)
    st_numscalar("mahal_zero",all(rdistnn_distance_matrix(J(3,2,7),(8,-1),"mahalanobis",MG):==0))
    st_numscalar("mahal_legacy",mreldif(rdistnn_distance_matrix(MB,(0,0),"mahalanobis",diag((.25,1))),ME))
}
test_assert "tiny covariance range" scalar(mahal_extreme_0)<2e-14
test_assert "huge covariance range" scalar(mahal_extreme_1)<2e-14
test_assert "rank-deficient Moore-Penrose geometry" scalar(mahal_singular)==1
test_assert "zero rank projects to zero distance" scalar(mahal_zero)==1
test_assert "legacy inverse geometry remains readable" scalar(mahal_legacy)<2e-15
mata {
    MX=(J(3,1,2^600),(-1\0\1)*2^-600)
    MG=rdistnn_fit_geometry(MX,"mahalanobis","sd",J(1,0,.),1e-12)
    st_numscalar("mahal_offset",MG[1,2]==1 & mreldif(rdistnn_distance_matrix(MX,(2^600,0),"mahalanobis",MG),(1\0\1))<2e-14)
    MX=((-1\0\1)*2^-600,J(3,1,0))
    MG=rdistnn_fit_geometry(MX,"mahalanobis","sd",J(1,0,.),1e-12)
    st_numscalar("mahal_null",mreldif(rdistnn_distance_matrix(MX,(0,2^1000),"mahalanobis",MG),(1\0\1)))
    MX=((-1\0\1)*2^1022,J(3,1,0))
    MG=rdistnn_fit_geometry(MX,"mahalanobis","sd",J(1,0,.),1e-12)
    st_numscalar("mahal_rawdiff",mreldif(rdistnn_distance_matrix(MX,(2^1022,0),"mahalanobis",MG),(2\1\0)))
    MX=(-1\0\1)
    MG=rdistnn_fit_geometry(MX,"mahalanobis","sd",J(1,0,.),1e-12)
    st_numscalar("mahal_one",mreldif(rdistnn_distance_matrix(MX,0,"mahalanobis",MG),(1\0\1)))
    MX=(-2^-900\2^-900)
    MG=rdistnn_fit_geometry(MX,"mahalanobis","sd",J(1,0,.),1e-12)
    st_numscalar("mahal_factor_range",mreldif(rdistnn_distance_matrix(MX,2^123,"mahalanobis",MG)/(2^1022*sqrt(2)),J(2,1,1)))
}
test_assert "constant level cannot erase small variation" scalar(mahal_offset)==1
test_assert "extreme null coordinate does not overflow projection" scalar(mahal_null)<2e-14
test_assert "scaled projection recovers overflowing raw difference" scalar(mahal_rawdiff)<2e-14
test_assert "existing one-covariate geometry remains supported" scalar(mahal_one)<2e-14
test_assert "factor-aware projection recovers finite distance" scalar(mahal_factor_range)<2e-14
mata {
    MX=(J(6,1,2^50):+(0\1\0\-1\0\0), (0\0\2^-40\0\-2^-40\rdistnn_power2(-1074)))
    MG=rdistnn_fit_geometry(MX,"mahalanobis","sd",J(1,0,.),1e-12)
    MD=rdistnn_distance_matrix(MX,MX[6,],"mahalanobis",MG)
    st_numscalar("mahal_subnormal",MD[6]==0 & MD[1]>0)
}
test_assert "offset fallback keeps nonzero subnormal target separation" scalar(mahal_subnormal)==1
