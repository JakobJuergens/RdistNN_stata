// Shared covariance arithmetic and estimator integration, with analytic values.
do "tests/test_framework.do"
test_section "Covariance numeric range"
mata {
    ca=1.25*2^511
    CX=(-ca,ca/2,7 \ ca,-ca/2,7)
    CV=rdistnn_covariance(CX,.5,1)
    CE=(ca,-ca/2,0)'*(ca,-ca/2,0)
    st_numscalar("cov_large",mreldif(CV/ca^2,CE/ca^2))
    st_numscalar("cov_sym",mreldif(CV,CV'))
    CT=J(8,1,rdistnn_power2(-538)):*((-1\1\-1\1\-1\1\-1\1))
    st_numscalar("cov_tiny",rdistnn_covariance(CT,7/8,1)[1,1]==rdistnn_power2(-1073))
    CX=((-2^-540,3,-1) \ (2^-540,3,1))
    CV=rdistnn_covariance(CX,.5,1,1)
    st_numscalar("cov_failure",all(CV[1,]:>=.) & all(CV[,1]:>=.) & CV[2,2]==0 & CV[3,3]==1)
    st_numscalar("cov_constant",rdistnn_covariance(J(3,1,2^900),1,1)[1,1]==0)
    CS=(2^1022,2^1021\2^1021,2^1022)
    st_numscalar("cov_midpoint",all(rdistnn_cov_symmetrize(CS):==CS))
}
test_assert "finite covariance after overflowing raw sum" scalar(cov_large)<2e-15
test_assert "exact symmetric covariance" scalar(cov_sym)==0
test_assert "aggregate subnormal variance recovered" scalar(cov_tiny)==1
test_assert "unavailable column preserves other covariance" scalar(cov_failure)==1
test_assert "large constant is deterministic" scalar(cov_constant)==1
test_assert "symmetrization preserves finite large diagonal" scalar(cov_midpoint)==1
mata {
    ca=rdistnn_power2(510); ct=rdistnn_power2(-510)
    CX=(ca,0 \ -ca,0 \ ct,ct \ -ct,-ct \ 0,ca \ 0,-ca)
    CV=rdistnn_covariance(CX,1,1)
    st_numscalar("cov_small_product",CV[1,2]/rdistnn_power2(-1019))
    ct=rdistnn_power2(-1000); cb=rdistnn_power2(500)
    CX=(ca,0 \ -ca,0 \ ct,cb \ -ct,-cb \ 0,ca \ 0,-ca)
    CV=rdistnn_covariance(CX,.5,1)
    st_numscalar("cov_small_input",CV[1,2]/rdistnn_power2(-500))
    CX=(ca,0 \ -ca,0 \ ct,cb \ 0,-cb \ 0,ca \ 0,-ca)
    CV=rdistnn_covariance(CX,5/6,1)
    st_numscalar("cov_small_center",CV[1,2]/rdistnn_power2(-500))
    st_numscalar("cov_small_permute",all(rdistnn_covariance(CX[,(2,1)],5/6,1):==CV[(2,1),(2,1)]))
}
test_assert "finite cross-products survive column scaling" abs(scalar(cov_small_product)-1)<3e-15
test_assert "finite original entries survive column scaling" abs(scalar(cov_small_input)-1)<3e-15
test_assert "original-scale centering preserves residual cross-products" abs(scalar(cov_small_center)-5/6)<3e-15
test_assert "exceptional cross-covariance respects column permutation" scalar(cov_small_permute)==1
mata {
    ca=1.5*rdistnn_power2(512)
    CX=(ca,ca \ -ca,ca \ ca,-ca \ -ca,-ca)
    CV=rdistnn_covariance(CX,1/64,1)
    st_numscalar("cov_protected",CV[1,2]==0 & max(abs(diagonal(CV)/(ca/4)^2:-1))<2e-15)
    ct=rdistnn_power2(-538)
    CX=(-ct,-ct \ ct,ct \ -ct,-ct \ ct,ct \ -ct,-ct \ ct,ct \ -ct,-ct \ ct,ct)
    st_numscalar("cov_underflow_protected",rdistnn_covariance(CX,7/8,1)[1,2]==rdistnn_power2(-1073))
    ca=rdistnn_power2(510); ct=rdistnn_power2(-1000); cb=rdistnn_power2(500)
    CX=(ca,0 \ -ca,0 \ ct,cb \ -ct,cb \ 0,ca \ 0,-ca)
    st_numscalar("cov_exact_zero",rdistnn_covariance(CX,1,1)[1,2]==0)
    CU=(-2,3 \ -1,-1 \ 0,2 \ 1,0 \ 2,-2 \ 3,1)
    CE=rdistnn_covariance(CU,5/6,1)
    CX=rdistnn_power2(511):+CU*rdistnn_power2(459)
    CV=rdistnn_covariance(CX,5/6,1)/rdistnn_power2(918)
    st_numscalar("cov_offset",mreldif(CV,CE))
}
test_assert "overflowing original products retain normalized recovery" scalar(cov_protected)==1
test_assert "underflowing original products retain normalized recovery" scalar(cov_underflow_protected)==1
test_assert "original-scale covariance retains exact cancellation" scalar(cov_exact_zero)==1
test_assert "exceptional covariance retains large common offsets" scalar(cov_offset)<3e-15
mata {
    ct=rdistnn_power2(-524); ca=1.1*ct; cb=ct*(1-rdistnn_power2(-27))
    CX=(ca,ct \ ca,-cb)
    CV=rdistnn_covariance(CX,1,0)
    st_numscalar("cov_subnormal_uncentered",CV[1,2]==rdistnn_power2(-1074))
    CX=(ca,ct \ ca,-cb \ -ca,-ct \ -ca,cb)
    CV=rdistnn_covariance(CX,1,1)
    st_numscalar("cov_subnormal_centered",CV[1,2]==rdistnn_power2(-1074))
    st_numscalar("cov_smallest",rdistnn_power2(-1074))
}
test_assert "normalized products retain a positive subnormal cross-entry" scalar(cov_subnormal_uncentered)==1
test_assert "normalized centered products retain a positive subnormal cross-entry" scalar(cov_subnormal_centered)==1
// Restore the genuine underflow fixture before testing its unavailable variance.
mata: CX=((-rdistnn_power2(-540),3,-1) \ (rdistnn_power2(-540),3,1))
capture mata: rdistnn_covariance(CX,.5,1)
test_assert "underflow rejects misleading zero uncertainty" _rc==498
foreach engine in mata plugin naive {
    global rdistnn_jackknife_engine `engine'
    mata {
        ca=1.25*2^511
        CV=rdistnn_jk_cov_multi_from_ranks((1\2),(-ca\ca),1,1,"exact")
        st_numscalar("cov_engine",CV[1,1]/ca^2)
    }
    test_assert "`engine' deletion finite covariance" abs(scalar(cov_engine)-1)<2e-15
    clear
    set obs 2
    gen double x=_n-1
    gen double y=cond(_n==1,-1,1)*1.25*2^511
    np_dnn y, xcols(x) x0(0) s(1)
    matrix CV=e(estimate_vcov)
    test_assert "`engine' NP finite covariance" abs(CV[1,1]/(1.25*2^511)^2-1)<2e-15
    replace y=cond(_n==1,-1,1)*2^-540
    capture np_dnn y, xcols(x) x0(0) s(1)
    test_assert "`engine' NP underflow error" _rc==498
    np_dnn y, xcols(x) x0(0) s(1) nose
    test_assert "`engine' point fit remains available" e(N)==2
    clear
    set obs 6
    gen double x=_n
    gen double numerator=cond(_n==1,2^513,cond(_n==2,-2^513,cond(_n==3,2^-507,cond(_n==4,-2^-507,0))))
    gen double denominator=cond(_n==5,2^513,cond(_n==6,-2^513,cond(_n==3,2^-507,cond(_n==4,-2^-507,0))))
    clate_dnn numerator, denominator(denominator) pseudooutcome xcols(x) x0(0) s(1)
    matrix CV=e(component_vcov)
    test_assert "`engine' paired component retains finite cross-covariance" abs(CV[1,2]/(2^-1014/15)-1)<2e-14
    clear
    set obs 4
    gen double x=_n
    gen double numerator=.
    gen double denominator=.
    mata {
        ct=rdistnn_power2(-524); ca=1.1*ct; cb=ct*(1-rdistnn_power2(-27))
        st_store(. ,"numerator",3*(ca\ca\-ca\-ca))
        st_store(. ,"denominator",3*(ct\-cb\-ct\cb))
    }
    clate_dnn numerator, denominator(denominator) pseudooutcome xcols(x) x0(0) s(1)
    matrix CV=e(component_vcov)
    test_assert "`engine' paired component retains normalized subnormal cross-entry" CV[1,2]==scalar(cov_smallest)
}
macro drop rdistnn_jackknife_engine
