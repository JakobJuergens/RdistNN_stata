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
}
macro drop rdistnn_jackknife_engine
