// End-to-end parity, including a default whose earlier duplicate precedes
// another computed order. Catch covariance reuse from the wrong candidate.
version 17
quietly do "src/_rdistnn_mata.ado"
do "tests/test_framework.do"
test_section "jackknife engines: native matrix bridge and fallbacks"
matrix JP_R=(1\2\3\4)
matrix JP_S=(2,1 \ 1,-1 \ 4,2 \ 3,1)
matrix JP_W=(2/3\1/3\0)
mata {
    JP_PAYLOAD=(1\4\1\2\0)\vec(st_matrix("JP_R"))\vec(st_matrix("JP_S"))\vec(st_matrix("JP_W"))
    st_matrix("JP_F",rdistnn_native_frame(JP_PAYLOAD))
    st_matrix("JP_D",rdistnn_native_frame(J(8,1,73)))
}
capture plugin call rdistnn_plugin, loo_packed JP_F JP_D
local native_rc=_rc
test_assert "native packed deletion bridge available" `native_rc'==0
if `native_rc'==0 {
    mata: st_matrix("JP_L",rdistnn_native_unframe(st_matrix("JP_D"),4,2)); st_numscalar("_je_literal",max(abs(st_matrix("JP_L")-(2,0\8/3,4/3\5/3,1/3\5/3,1/3))))
    test_assert "native bridge preserves paired literal deletions" scalar(_je_literal)<1e-12
    foreach fault in duplicate noninteger missing {
        matrix JP_BAD=(1\2\3\4)
        if "`fault'"=="duplicate" matrix JP_BAD[2,1]=1
        if "`fault'"=="noninteger" matrix JP_BAD[4,1]=4.5
        if "`fault'"=="missing" matrix JP_BAD[4,1]=.
        mata: JP_PAYLOAD=(1\4\1\2\0)\vec(st_matrix("JP_BAD"))\vec(st_matrix("JP_S"))\vec(st_matrix("JP_W")); st_matrix("JP_F",rdistnn_native_frame(JP_PAYLOAD)); st_matrix("JP_D",rdistnn_native_frame(J(8,1,73))); JP_D_BEFORE=st_matrix("JP_D")
        capture plugin call rdistnn_plugin, loo_packed JP_F JP_D
        local bad_rc=_rc
        mata: st_numscalar("_je_unchanged",all(st_matrix("JP_D"):==JP_D_BEFORE))
        test_assert "native bridge rejects `fault' ranks before output writes" `bad_rc'==198 & scalar(_je_unchanged)
    }
    capture plugin call rdistnn_plugin, loo_packed JP_F
    test_assert "native bridge rejects missing matrix argument" _rc==198
}

global rdistnn_jackknife_engine auto
global rdistnn_use_plugin 0
mata: _je_mata=rdistnn_loo_from_ranks(st_matrix("JP_R"),st_matrix("JP_S"),1,2,"exact")
global rdistnn_use_plugin 1
global rdistnn_plugin_ok 0
capture mata: _je_L=rdistnn_loo_from_ranks(st_matrix("JP_R"),st_matrix("JP_S"),1,2,"exact")
test_assert "required unavailable plugin fails explicitly" _rc!=0
capture macro drop rdistnn_use_plugin
mata: st_numscalar("_je_fallback",max(abs(rdistnn_loo_from_ranks(st_matrix("JP_R"),st_matrix("JP_S"),1,2,"exact")-_je_mata)))
test_assert "unavailable plugin falls back to fast Mata" scalar(_je_fallback)<1e-12
macro drop rdistnn_plugin_ok
quietly do "src/_rdistnn_mata.ado"

// Mata-created bridge matrices can exceed the user's Stata matsize setting.
// Exercise both tall inputs and output widths exceeding it without changing it.
local bridge_matsize=c(matsize)
mata {
    _je_n=st_numscalar("c(matsize)")+1
    _je_r=(1::_je_n)
    _je_s=sin(_je_r)
    _je_wideR=(1::4)*J(1,ceil(_je_n/2),1)
    _je_wideS=(2,1 \ 1,-1 \ 4,2 \ 3,1)
}
foreach method in exact asymp {
    global rdistnn_jackknife_engine mata
    mata: _je_tallref=rdistnn_loo_from_ranks(_je_r,_je_s,1,2,"`method'"); _je_wideref=rdistnn_loo_from_ranks(_je_wideR,_je_wideS,1,2,"`method'")
    foreach engine in plugin auto {
        global rdistnn_jackknife_engine `engine'
        global rdistnn_use_plugin 1
        foreach shape in tall wide {
            if "`shape'"=="tall" local inputs "_je_r,_je_s"
            else local inputs "_je_wideR,_je_wideS"
            capture mata: _je_native=rdistnn_loo_from_ranks(`inputs',1,2,"`method'")
            local bridge_rc=_rc
            test_assert "`engine'/`method': native `shape' matrices above matsize succeed" `bridge_rc'==0
            if `bridge_rc'==0 {
                mata: st_numscalar("_je_large",max(abs(_je_native-_je_`shape'ref)))
                test_assert "`engine'/`method': native `shape' deletions above matsize match Mata" scalar(_je_large)<1e-12
            }
        }
    }
}
test_assert "native bridge preserves user's matsize setting" c(matsize)==`bridge_matsize'
macro drop rdistnn_use_plugin

test_section "jackknife engines: estimator covariance and candidate reuse"
clear
set obs 24
gen double x=(_n-12)/5
gen double y=sin(x)+x/3
gen double den=1+cos(x)/5
matrix JE_X0=(-.5\0\.7)
foreach method in exact asymp {
    foreach estimator in np cate clate {
        local extra ""
        if "`estimator'"=="cate" local extra "pseudooutcome"
        if "`estimator'"=="clate" local extra "pseudooutcome denominator(den)"
        global rdistnn_jackknife_engine naive
        quietly `estimator'_dnn y, xcols(x) x0matrix(JE_X0) k(2) ///
            method(`method') tuningmode(grid) lambdagrid(.37 .75 .4 .5) `extra'
        matrix JE_B=e(b)
        matrix JE_V=e(V)
        matrix JE_G=e(tuning_se)
        test_assert "`estimator'/`method': repeated default follows another order" e(tuning_default_candidate)==3
        foreach engine in mata plugin auto {
            if "`engine'"=="plugin" & `native_rc'!=0 continue
            global rdistnn_jackknife_engine `engine'
            quietly `estimator'_dnn y, xcols(x) x0matrix(JE_X0) k(2) ///
                method(`method') tuningmode(grid) lambdagrid(.37 .75 .4 .5) `extra'
            mata: st_numscalar("_je_err",max((max(abs(st_matrix("JE_B")-st_matrix("e(b)"))), max(abs(st_matrix("JE_V")-st_matrix("e(V)"))), max(abs(st_matrix("JE_G")-st_matrix("e(tuning_se)"))))))
            test_assert "`estimator'/`method'/`engine': fit, joint covariance and grid SEs match naive" scalar(_je_err)<1e-10
        }
    }
}

rdnn_dgp_plqr, n(40) seed(8401)
foreach method in exact asymp {
    global rdistnn_jackknife_engine naive
    quietly plqr_dnn y, x(x1 x2) x0(0 0) d(d1 d2) g(q) md(md1 md2) ///
        thetaeval(1 -.5) tau(.5) k(2) s(6) method(`method')
    matrix JE_V=e(moment_vcov)
    foreach engine in mata plugin auto {
        if "`engine'"=="plugin" & `native_rc'!=0 continue
        global rdistnn_jackknife_engine `engine'
        quietly plqr_dnn y, x(x1 x2) x0(0 0) d(d1 d2) g(q) md(md1 md2) ///
            thetaeval(1 -.5) tau(.5) k(2) s(6) method(`method')
        mata: st_numscalar("_je_err",max(abs(st_matrix("JE_V")-st_matrix("e(moment_vcov)"))))
        test_assert "PLQR/`method'/`engine': vector covariance matches naive" scalar(_je_err)<1e-10
    }
}
macro drop rdistnn_jackknife_engine
