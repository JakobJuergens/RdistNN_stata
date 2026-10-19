version 17
do "tests/test_framework.do"
quietly do "src/_rdistnn_mata.ado"
test_section "translation-safe numeric centering"
foreach engine in mata plugin naive {
    global rdistnn_jackknife_engine `engine'
    clear
    set obs 50
    gen double fx=_n
    gen double fy=1e16+2*mod(_n+1,2)
    replace fy=-1e16 in 50
    gen double fc=fy-1e16
    gen double fd=1
    foreach family in np cate clate {
        local extra ""
        if "`family'"=="cate" local extra "pseudooutcome"
        if "`family'"=="clate" local extra "pseudooutcome denominator(fd)"
        foreach outcome in fy fc {
            quietly `family'_dnn `outcome', xcols(fx) x0(0) s(49) method(asymp) `extra'
            matrix FV=e(V)
            test_assert "`engine' `family' `outcome': exact rational SE anchor" abs(sqrt(FV[1,1])-.983985898337929651)<2e-14
        }
    }
}
macro drop rdistnn_jackknife_engine
mata:
_fn_sd=rdistnn_sample_sd(1e16:+(0\2\0\2))
_fn_cov=rdistnn_covariance((1e16:+(0\2\0\2),7:+3*(0\2\0\2)),1/3,1)
st_numscalar("_fn_sd",abs(_fn_sd-sqrt(4/3))<1e-14)
st_numscalar("_fn_cov",max(abs(_fn_cov-(4/3,4\4,12)))<1e-13)
end
test_assert "observed-origin SD preserves ordinary local variation" scalar(_fn_sd)
test_assert "observed-origin covariance preserves ordinary local variation" scalar(_fn_cov)
foreach engine in mata plugin {
    global rdistnn_distance_engine `engine'
    clear
    set obs 4
    gen double fx=cond(mod(_n,2),0,2)
    gen double fz=(_n>2)
    gen double fy=10*_n
    quietly np_dnn fy, xcols(fx fz) x0(0 .9) s(3) distance(l2_scaled) nose
    matrix F0=e(b)
    replace fx=fx+1e16
    quietly np_dnn fy, xcols(fx fz) x0(1e16 .9) s(3) distance(l2_scaled) nose
    matrix F1=e(b)
    test_assert "`engine': translated SD geometry preserves estimate" F0[1,1]==F1[1,1]
}
macro drop rdistnn_distance_engine
// Input row permutation must not move the numerical origin to a remote tail.
foreach engine in mata plugin naive {
    global rdistnn_jackknife_engine `engine'
    clear
    set obs 50
    gen double fx=51-_n
    gen double fy=1e16+2*mod(fx+1,2)
    replace fy=-1e16 if fx==50
    quietly np_dnn fy, xcols(fx) x0(0) s(49) method(asymp)
    matrix FV=e(V)
    test_assert "`engine': reordered remote tail preserves local SE" abs(sqrt(FV[1,1])-.983985898337929651)<2e-14
}
macro drop rdistnn_jackknife_engine
mata: _fn_C=rdistnn_covariance((2^55,1\3,-1\-2^55,1\17,-1),1,1)
mata: st_numscalar("_fn_cross",_fn_C[1,2])
test_assert "covariance origin does not erase represented small cross terms" scalar(_fn_cross)==-20
mata {
    _fn_Z=(2^55-4,1 \ -2,-1 \ -2^55-8,1 \ 12,-1)
    _fn_L=(-1,0 \ 0,0 \ 3,0 \ 0,0)
    _fn_C=rdistnn_cov_crossprod(_fn_Z,_fn_L)
    st_numscalar("_fn_residual",_fn_C[1,2])
    _fn_C=rdistnn_cov_crossprod((2^55,1 \ 3,-1 \ -2^55,1 \ 17,-1),J(0,0,0))
    st_numscalar("_fn_raw",_fn_C[1,2])
    // Mata's ^ can round these scales; use exact binary inputs for this identity.
    _fn_C=rdistnn_covariance((2^55,1 \ 3,-1 \ -2^55,1 \ 17,-1):*(rdistnn_power2(454),rdistnn_power2(-450)),1,1)
    st_numscalar("_fn_scaled",_fn_C[1,2])
}
test_assert "sensitive covariance retains mean-subtraction residuals" scalar(_fn_residual)==-20
test_assert "uncentered cross-products retain small represented terms" scalar(_fn_raw)==-20
test_assert "sensitive covariance retains residuals after binary scaling" scalar(_fn_scaled)==-320
