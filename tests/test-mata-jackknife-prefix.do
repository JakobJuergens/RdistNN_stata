// All engines must retain observation identity and paired/cross-target signals.
version 17
quietly do "src/_rdistnn_mata.ado"
do "tests/test_framework.do"

test_section "prefix/suffix jackknife: literal paired deletions"
global rdistnn_jackknife_engine mata
capture mata: _jp_L = rdistnn_loo_from_ranks((1::4), (2,1 \ 1,-1 \ 4,2 \ 3,1), 1, 2, "exact")
local prefix_rc = _rc
test_assert "Mata prefix engine is available" `prefix_rc' == 0
if `prefix_rc' == 0 {
    mata: st_numscalar("_jp_literal", max(abs(_jp_L-(2,0 \ 8/3,4/3 \ 5/3,1/3 \ 5/3,1/3))))
    test_assert "paired deletions match hand calculation" scalar(_jp_literal) < 1e-12
}

// Explicit deletion oracle, independently retaining the preexisting algorithm.
capture mata: mata drop jp_explicit()
mata:
real matrix jp_explicit(real matrix ranks, real matrix signal,
    real scalar k, real scalar s, string scalar method)
{
    real scalar n, m, q, i, j
    real colvector keep, r, w
    real matrix out
    n=rows(ranks); m=cols(ranks); q=cols(signal)
    out=J(n,m*q,.)
    for(i=1;i<=n;i++) {
        keep=((1::n):!=i)
        for(j=1;j<=m;j++) {
            r=select(ranks[,j],keep)
            r=r:-(r:>ranks[i,j])
            w=rdistnn_weights_from_ranks(r,k,s,method)
            out[i,((j-1)*q+1)..(j*q)]=cross(w,select(signal,keep))
        }
    }
    return(out)
}
end

if `prefix_rc' == 0 {
    foreach method in exact asymp {
        foreach engine in mata naive auto plugin {
            global rdistnn_jackknife_engine `engine'
            mata {
                _jp_R=(1,6 \ 3,4 \ 5,2 \ 6,1 \ 4,3 \ 2,5)
                _jp_S=(1,1,7 \ 3,-2,7 \ -2,1,7 \ 4,3,7 \ 2,2,7 \ 5,1,7)
                _jp_err=0
                for(_jp_k=1;_jp_k<=3;_jp_k++) {
                    for(_jp_s=_jp_k;_jp_s<=5;_jp_s++) {
                        _jp_L=rdistnn_loo_from_ranks(_jp_R,_jp_S,_jp_k,_jp_s,"`method'")
                        _jp_want=jp_explicit(_jp_R,_jp_S,_jp_k,_jp_s,"`method'")
                        _jp_err=max((_jp_err,max(abs(_jp_L-_jp_want))))
                    }
                }
                st_numscalar("_jp_error",_jp_err)
            }
            test_assert "`engine'/`method': every deletion, component and target matches explicit reference" scalar(_jp_error)<1e-11
        }
    }
    global rdistnn_jackknife_engine mata
    mata: st_numscalar("_jp_n2",max(abs(rdistnn_loo_from_ranks((2\1),(9\4),1,1,"exact")-(4\9))))
    test_assert "two-observation endpoint" scalar(_jp_n2)<1e-12
    capture mata: _jp_L=rdistnn_loo_from_ranks((1\1),(2\3),1,1,"exact")
    test_assert "duplicate ranks rejected" _rc!=0
    capture mata: _jp_L=rdistnn_loo_from_ranks((1::4),(1::4),1,4,"exact")
    test_assert "s=n rejected for deletions" _rc!=0

    // An unavailable fast Mata engine must retain a functional final fallback.
    global rdistnn_jackknife_engine auto
    global rdistnn_use_plugin 0
    mata: mata drop rdistnn_loo_prefix()
    mata: st_numscalar("_jp_last",max(abs(rdistnn_loo_from_ranks((1::4),(2\1\4\3),1,2,"exact")-(2\8/3\5/3\5/3))))
    test_assert "missing prefix engine falls back to naive deletions" scalar(_jp_last)<1e-12
    quietly do "src/_rdistnn_mata.ado"
}
macro drop rdistnn_jackknife_engine

test_section "jackknife centering: asymptotic tails and unused extremes"
foreach engine in naive mata plugin auto {
    global rdistnn_jackknife_engine `engine'
    mata {
        _jp_tail=(0\0\J(8,1,1))
        _jp_A=(8/9)^9
        _jp_B=(7/9)^9
        _jp_expected=(36/25)*(_jp_A-_jp_B)^2
        _jp_L=rdistnn_loo_from_ranks((1::10),_jp_tail,1,9,"asymp")
        _jp_V=rdistnn_jk_cov_from_ranks((1::10),_jp_tail,1,9,"asymp")
        st_numscalar("_jp_tail_error",max(abs(_jp_L-(J(2,1,_jp_A)\J(8,1,_jp_B)))))
        st_numscalar("_jp_tail_variance",abs(_jp_V-_jp_expected))
        st_numscalar("_jp_exact_zero",rdistnn_jk_cov_from_ranks((1::10),_jp_tail,1,9,"exact"))
    }
    test_assert "`engine': asymptotic deletions retain weight beyond the exact cutoff" scalar(_jp_tail_error)<1e-12
    test_assert "`engine': asymptotic variance matches independent R/Stata fixture" scalar(_jp_tail_variance)<1e-12
    test_assert "`engine': exact rule excludes the same tail" scalar(_jp_exact_zero)==0
    foreach sign in 1 -1 {
        // The last row's subtraction would exceed Stata's finite-value range.
        mata: st_numscalar("_jp_extreme",rdistnn_jk_cov_from_ranks((1::30),`sign'*(J(29,1,6e307)\-6e307),1,10,"exact"))
        test_assert "`engine'/sign=`sign': unused extreme cannot overflow centering" scalar(_jp_extreme)==0
    }
}
macro drop rdistnn_jackknife_engine

test_section "jackknife covariance: zero-weight opposite-sign observations"
foreach engine in naive mata plugin auto {
    global rdistnn_jackknife_engine `engine'
    foreach k in 1 3 {
        foreach sign in 1 -1 {
            mata {
                _jp_z=`sign'*(J(29,1,1e16)\-1)
                _jp_v=rdistnn_jk_cov_from_ranks((1::30),_jp_z,`k',10,"exact")
                st_numscalar("_jp_support",_jp_v)
                _jp_raw=rdistnn_loo_from_ranks((1::30),_jp_z,`k',10,"exact")
                st_numscalar("_jp_rawspread",max(_jp_raw)-min(_jp_raw))
            }
            test_assert "`engine'/k=`k'/sign=`sign': unused opposite-sign row cannot create variance" scalar(_jp_support)==0
            test_assert "`engine'/k=`k'/sign=`sign': raw deletions remain constant" scalar(_jp_rawspread)==0
        }
        mata {
            _jp_base=4*(1::30)
            _jp_shifted=1e16:+_jp_base
            _jp_shifted[30]=-1e307
            _jp_v0=rdistnn_jk_cov_from_ranks((1::30),_jp_base,`k',10,"exact")
            _jp_v1=rdistnn_jk_cov_from_ranks((1::30),_jp_shifted,`k',10,"exact")
            st_numscalar("_jp_local",abs(_jp_v1-_jp_v0))
        }
        test_assert "`engine'/k=`k': large unused row preserves local variation" scalar(_jp_local)<1e-9
    }
    mata {
        _jp_r=(1::60)
        _jp_R=(_jp_r,61:-_jp_r)
        _jp_z=(4*(1::30)\-8*(1::30))
        _jp_shifted=_jp_z+(J(30,1,1e16)\J(30,1,-1e16))
        _jp_V0=rdistnn_jk_cov_multi_from_ranks(_jp_R,_jp_z,3,40,"exact")
        _jp_V1=rdistnn_jk_cov_multi_from_ranks(_jp_R,_jp_shifted,3,40,"exact")
        st_numscalar("_jp_target",max(abs(_jp_V1-_jp_V0)))
    }
    test_assert "`engine': target-specific supports preserve joint covariance" scalar(_jp_target)<1e-9
    // Shifted deletions (4,0,0,0) give V=(3/4)*(9+1+1+1)=9;
    // this fixture is also checked by the R package.
    mata: st_numscalar("_jp_golden",rdistnn_jk_cov_from_ranks((1::4),(1e16:+(4\8\12))\-1,1,3,"exact"))
    test_assert "`engine': large-level R/Stata golden variance is nine" abs(scalar(_jp_golden)-9)<1e-12
}
macro drop rdistnn_jackknife_engine

test_section "jackknife support centering: all estimator commands"
clear
set obs 30
gen double jp_x=_n
gen double jp_y=cond(_n<30,1e16,-1)
gen double jp_den=1
gen double jp_d=2*jp_y
gen double jp_zero=0
gen double jp_positive=1
foreach engine in naive mata plugin auto {
    global rdistnn_jackknife_engine `engine'
    foreach estimator in np cate clate {
        local extra ""
        if "`estimator'"=="cate" local extra "pseudooutcome"
        if "`estimator'"=="clate" local extra "pseudooutcome denominator(jp_den)"
        quietly `estimator'_dnn jp_y, xcols(jp_x) x0(0) k(3) s(10) `extra'
        matrix JP_SUPPORT_V=e(V)
        test_assert "`engine'/`estimator': unused opposite-sign row leaves zero variance" JP_SUPPORT_V[1,1]==0
    }
    quietly plqr_dnn jp_positive, x(jp_x) x0(0) d(jp_d) g(jp_zero) md(jp_zero) thetaeval(0) tau(.5) k(3) s(10)
    matrix JP_SUPPORT_V=e(moment_vcov)
    test_assert "`engine'/PLQR: unused opposite-sign row leaves zero moment variance" JP_SUPPORT_V[1,1]==0
}
macro drop rdistnn_jackknife_engine
capture macro drop rdistnn_use_plugin

// Covariance is invariant to a common signal level. Raw split sums at 1e16
// can manufacture variation even for exactly constant outcomes.
test_section "jackknife covariance: constants and large common levels"
quietly do "src/_rdistnn_clate_mata.ado"
foreach engine in naive mata plugin auto {
    global rdistnn_jackknife_engine `engine'
    mata {
        _jp_i=(1::30)
        _jp_R=(_jp_i,31:-_jp_i)
        _jp_S=(J(30,1,1e16),1e16:+4*_jp_i)
        _jp_V=rdistnn_jk_cov_from_ranks(_jp_i,_jp_S,1,1,"exact")
        _jp_M=rdistnn_jk_cov_multi_from_ranks(_jp_R,_jp_S[,2],1,1,"exact")
        st_numscalar("_jp_const",max(abs(_jp_V[1,])))
        st_numscalar("_jp_shift",max((abs(_jp_V[2,2]-124/3),max(abs(_jp_M:-124/3)))))
    }
    test_assert "`engine': exactly constant large signal has zero covariance" scalar(_jp_const)==0
    test_assert "`engine': large signal level preserves sample-mean variance and joint covariance" scalar(_jp_shift)<1e-9
    mata {
        _jp_C=_jp_CV=_jp_TV=_jp_status=J(0,0,.)
        _jp_theta=J(1,0,.)
        rdistnn_clate_fit(_jp_R,(_jp_S[,2],J(30,1,1)),1,1,"exact",1,1e-8,
            _jp_C,_jp_theta,_jp_TV,_jp_CV,_jp_status)
        st_numscalar("_jp_clate_shift",max(abs(_jp_TV:-124/3)))
        rdistnn_clate_fit(_jp_R,(J(30,1,1e16),J(30,1,1)),1,1,"exact",1,1e-8,
            _jp_C,_jp_theta,_jp_TV,_jp_CV,_jp_status)
        st_numscalar("_jp_clate_const",max(abs(_jp_TV)))
    }
    test_assert "`engine': paired CLATE covariance preserves large-level variation" scalar(_jp_clate_shift)<1e-9
    test_assert "`engine': constant paired CLATE signals have zero variance" scalar(_jp_clate_const)==0
    foreach sign in 1 -1 {
        mata: st_numscalar("_jp_outlier",abs(rdistnn_jk_cov_from_ranks((1::4),`sign'*(1\2\3\1e20),1,3,"exact")-9/16))
        test_assert "`engine': irrelevant signed outlier preserves local variation (`sign')" scalar(_jp_outlier)<1e-12
    }
}
macro drop rdistnn_jackknife_engine
