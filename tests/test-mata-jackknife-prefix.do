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
