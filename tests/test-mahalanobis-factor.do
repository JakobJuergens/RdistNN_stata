// The true nearest neighbor is known in latent orthogonal coordinates.
do "tests/test_framework.do"
test_section "Factored Mahalanobis geometry"
test_section "Computed negative extreme retains finite fitted geometry"
clear
set obs 3
gen double mneg_x=_n*1e-300
replace mneg_x=(-5e307)*2 in 1
gen double mneg_zero=0
gen double mneg_y=_n
mata {
    MNEG=st_data(.,"mneg_x mneg_zero")
    MNEG_T=(0,0 \ MNEG[1,])
    st_matrix("MNEG_T",MNEG_T)
    st_numscalar("mneg_sd",rdistnn_sample_sd(MNEG[,1])/(-(MNEG[1,1]/sqrt(3))))
}
test_assert "computed negative sample SD matches independent rank-one reference" abs(scalar(mneg_sd)-1)<2e-14
foreach engine in mata plugin {
    global rdistnn_distance_engine `engine'
    foreach metric in l2_scaled mahalanobis {
        mata {
            MNEG_G=rdistnn_fit_geometry(MNEG,"`metric'","sd",J(1,0,.),1e-12)
            MNEG_D=rdistnn_distance_matrix(MNEG,MNEG_T,"`metric'",MNEG_G)
            MNEG_E=(sqrt(3),0 \ 0,sqrt(3) \ 0,sqrt(3))
            st_numscalar("mneg_distance",mreldif(MNEG_D,MNEG_E))
        }
        test_assert "`engine' `metric': negative finite extreme has analytic distances" scalar(mneg_distance)<2e-14
        if "`metric'"=="mahalanobis" {
            mata: st_numscalar("mneg_rank",MNEG_G[1,2])
            test_assert "`engine': negative extreme retains rank one" scalar(mneg_rank)==1
        }
        quietly weights_dnn mneg_x mneg_zero, x0matrix(MNEG_T) k(1) s(3) distance(`metric')
        matrix MNEG_W=r(weights)
        matrix MNEG_R=r(ranks)
        test_assert "`engine' `metric': exact independent nearest ranks" MNEG_R[1,1]==3 & MNEG_R[2,1]==1 & MNEG_R[3,1]==2 & MNEG_R[1,2]==1 & MNEG_R[2,2]==2 & MNEG_R[3,2]==3
        test_assert "`engine' `metric': requested engine executed" "$rdistnn_last_distance_engine"=="`engine'"
        quietly np_dnn mneg_y, xcols(mneg_x mneg_zero) x0matrix(MNEG_T) k(1) s(3) distance(`metric') nose
        matrix MNEG_F=e(estimates)
        test_assert "`engine' `metric': supported public NP averages both exact zero-distance neighbors" MNEG_F[1,1]==2.5 & MNEG_F[1,2]==1
        tempfile mnegsaved
        estimates save `mnegsaved', replace
        preserve
        clear
        estimates use `mnegsaved'
        quietly rdistnn_predict, x0matrix(MNEG_T)
        matrix MNEG_P=r(fit)
        test_assert "`engine' `metric': current saved prediction retains finite geometry" MNEG_P[1,1]==2.5 & MNEG_P[1,2]==1
        restore
        mata {
            mneg_at=(3\1\2)
            MNEG_PG=rdistnn_fit_geometry(MNEG[mneg_at,],"`metric'","sd",J(1,0,.),1e-12)
            MNEG_PD=rdistnn_distance_matrix(MNEG[mneg_at,],MNEG_T,"`metric'",MNEG_PG)
            st_numscalar("mneg_permutation",mreldif(MNEG_PD,MNEG_E[mneg_at,]))
        }
        test_assert "`engine' `metric': negative geometry survives row permutation" scalar(mneg_permutation)<2e-14
    }
}
capture macro drop rdistnn_distance_engine
clear
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
    MOX=(rdistnn_power2(100)\1\2\3)
    MOT=(2.1\1.1\2.1)
    MOE=J(4,3,.)
    for(mj=1;mj<=3;mj++) MOE[,mj]=abs(MOX:-MOT[mj])/rdistnn_power2(99)
    MOR=(4,4,4 \ 3,1,3 \ 1,2,1 \ 2,3,2)
    MOG=rdistnn_fit_geometry(MOX,"mahalanobis","sd",J(1,0,.),1e-12)
    st_matrix("MOX",MOX); st_matrix("MOT",MOT)
}
foreach engine in mata plugin {
    global rdistnn_distance_engine `engine'
    mata {
        MOD=rdistnn_distance_matrix(MOX,MOT,"mahalanobis",MOG)
        st_numscalar("mahal_origin_distances",max(abs(MOD:/MOE:-1)))
        st_numscalar("mahal_origin_ranks",all(rdistnn_ranks_from_dmat(MOD):==MOR))
        mshape=1
        for(mj=1;mj<=3;mj++) mshape=mshape & all(rdistnn_distance_matrix(MOX,MOT[mj],"mahalanobis",MOG):==MOD[,mj])
        st_numscalar("mahal_origin_shape",mshape)
        moi=(2\3\4\1)
        MOP=rdistnn_fit_geometry(MOX[moi],"mahalanobis","sd",J(1,0,.),1e-12)
        MOPD=rdistnn_distance_matrix(MOX[moi],MOT,"mahalanobis",MOP)
        st_numscalar("mahal_origin_permute",max(abs(MOPD:/MOE[moi,]:-1)))
    }
    test_assert "`engine': small nonzero Mahalanobis distances" scalar(mahal_origin_distances)<1e-14
    test_assert "`engine': independent exact small-neighbor ranks" scalar(mahal_origin_ranks)==1
    test_assert "`engine': small distances retain scalar/batch identity" scalar(mahal_origin_shape)==1
    test_assert "`engine': small distances survive training-row permutation" scalar(mahal_origin_permute)<1e-14
    test_assert "`engine': small-neighbor distance engine ran" "$rdistnn_last_distance_engine"=="`engine'"
    clear
    svmat double MOX, names(mox)
    gen double moy=cond(_n==1,100,(_n-1)*10)
    quietly np_dnn moy, xcols(mox1) x0matrix(MOT) k(1) s(4) distance(mahalanobis) nose
    matrix MOFIT=e(estimates)
    test_assert "`engine': exact nearest-neighbor NP estimates" MOFIT[1,1]==20 & MOFIT[1,2]==10 & MOFIT[1,3]==20
    tempfile mosaved
    estimates save `mosaved', replace
    clear
    estimates use `mosaved'
    matrix MONEW=MOT[2,1...] \ MOT[1,1...] \ MOT[2,1...]
    quietly rdistnn_predict, x0matrix(MONEW)
    matrix MOFIT=r(fit)
    test_assert "`engine': frozen small-neighbor saved prediction" MOFIT[1,1]==10 & MOFIT[1,2]==20 & MOFIT[1,3]==10
}
macro drop rdistnn_distance_engine
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
}
test_assert "tiny covariance range" scalar(mahal_extreme_0)<2e-14
test_assert "huge covariance range" scalar(mahal_extreme_1)<2e-14
test_assert "rank-deficient Moore-Penrose geometry" scalar(mahal_singular)==1
test_assert "zero rank projects to zero distance" scalar(mahal_zero)==1
capture mata: rdistnn_distance_matrix(MB,(0,0),"mahalanobis",diag((.25,1)))
test_assert "historical inverse geometry requires a fresh fit" _rc==498
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
