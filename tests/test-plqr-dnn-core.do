// test-plqr-dnn-core.do
// Tests for plqr_dnn: smoke behavior, eclass results, and core identities.
//
// Tags: [plqr] [smoke] [eclass] [identity]

quietly do "src/_rdistnn_mata.ado"
do "tests/test_framework.do"

// ════════════════════════════════════════════════════════════
test_section "plqr_dnn: [smoke] basic run and eclass results"
// ════════════════════════════════════════════════════════════

rdnn_dgp_plqr, n(120) seed(401)
plqr_dnn y, x(x1 x2) x0(0 0) d(d1 d2) g(q) md(md1 md2) ///
    thetaeval(1 -0.5) tau(0.5) k(4) s(12)

matrix B = e(b)
matrix V = e(V)
matrix Moment = e(moment)
matrix MomentV = e(moment_vcov)
matrix ThetaEval = e(theta_eval)
local cmd "`e(cmd)'"
local theta_source "`e(theta_source)'"
local inference_target "`e(inference_target)'"
local root_solver "`e(root_solver)'"
local kd = e(kd)
local tau = e(tau)
local nb = colsof(B)
local bcols : colnames B
local brows : rownames V
local bcolsV : colnames V

scalar _b1 = B[1,1]
scalar _b2 = B[1,2]
scalar _v11 = V[1,1]
scalar _v22 = V[2,2]
scalar _v12 = V[1,2]
scalar _v21 = V[2,1]
scalar _theta1 = ThetaEval[1,1]
scalar _theta2 = ThetaEval[1,2]

test_assert "[smoke] e(cmd) = plqr_dnn" "`cmd'" == "plqr_dnn"
mata: st_numscalar("_moment_alias", max(abs(st_matrix("B") - st_matrix("Moment"))))
mata: st_numscalar("_vcov_alias", max(abs(st_matrix("V") - st_matrix("MomentV"))))
test_assert "[smoke] e(moment) is the explicit moment copy of e(b)" scalar(_moment_alias) < 1e-12
test_assert "[smoke] e(moment_vcov) is the explicit covariance copy of e(V)" scalar(_vcov_alias) < 1e-12
test_assert "[smoke] metadata marks supplied theta and moment inference" ///
    "`theta_source'" == "supplied" & "`inference_target'" == "moment" & ///
    "`root_solver'" == "none"
test_assert "[smoke] e(b) has two moment components" `nb' == 2
test_assert "[smoke] e(b) column names follow d()" "`bcols'" == "d1 d2"
test_assert "[smoke] e(V) row names follow d()" "`brows'" == "d1 d2"
test_assert "[smoke] e(V) column names follow d()" "`bcolsV'" == "d1 d2"
test_assert "[smoke] moments are finite" scalar(_b1) < . & scalar(_b2) < .
test_assert "[smoke] variances are finite and non-negative" ///
    scalar(_v11) < . & scalar(_v11) >= 0 & scalar(_v22) < . & scalar(_v22) >= 0
test_assert "[smoke] covariance matrix is symmetric" abs(scalar(_v12) - scalar(_v21)) < 1e-10
test_assert "[smoke] e(theta_eval) stores supplied theta" ///
    abs(scalar(_theta1) - 1) < 1e-12 & abs(scalar(_theta2) + 0.5) < 1e-12
test_assert "[smoke] e(kd) counts d() variables" `kd' == 2
test_assert "[smoke] e(tau) stores tau" abs(`tau' - 0.5) < 1e-12

cap noisily plqr_dnn
local rc_replay = _rc
test_assert "[smoke] replay runs after plqr_dnn" `rc_replay' == 0

// ════════════════════════════════════════════════════════════
test_section "plqr_dnn: [identity] weighted moment equals manual sum"
// ════════════════════════════════════════════════════════════

rdnn_dgp_plqr, n(110) seed(402)
weights_dnn x1 x2, x0(0 0) k(4) s(10) generate(w_plqr)
gen double g1 = (0.5 - (y <= d1 - 0.5 * d2 + q)) * (d1 - md1)
gen double g2 = (0.5 - (y <= d1 - 0.5 * d2 + q)) * (d2 - md2)
gen double wg1 = w_plqr * g1
gen double wg2 = w_plqr * g2
quietly summarize wg1, meanonly
local m1 = r(sum)
quietly summarize wg2, meanonly
local m2 = r(sum)

plqr_dnn y, x(x1 x2) x0(0 0) d(d1 d2) g(q) md(md1 md2) ///
    thetaeval(1 -0.5) tau(0.5) k(4) s(10)
matrix B = e(moment)
scalar _r1 = abs(B[1,1] - `m1')
scalar _r2 = abs(B[1,2] - `m2')
test_assert "[identity] first moment equals manual weighted sum" scalar(_r1) < 1e-8
test_assert "[identity] second moment equals manual weighted sum" scalar(_r2) < 1e-8

// ════════════════════════════════════════════════════════════
test_section "plqr_dnn: [identity] covariance equals direct Mata helper"
// ════════════════════════════════════════════════════════════

rdnn_dgp_plqr, n(80) seed(403)
plqr_dnn y, x(x1 x2) x0(0.2 -0.1) d(d1 d2) g(q) md(md1 md2) ///
    thetaeval(1 -0.5) tau(0.5) k(3) s(8)
matrix V_cmd = e(moment_vcov)

mata {
    _X     = st_data(., "x1 x2")
    _D     = st_data(., "d1 d2")
    _md    = st_data(., "md1 md2")
    _q     = st_data(., "q")
    _y     = st_data(., "y")
    _beta  = (1 \ -0.5)
    _x0    = (0.2, -0.1)
    _ranks = rdistnn_localization_ranks(rdistnn_distance_matrix(_X,_x0,"l2"))
    _resid = 0.5 :- (_y :<= (_D * _beta :+ _q))
    _g     = (_D :- _md) :* (_resid * J(1, cols(_D), 1))
    _V_dir = rdistnn_jk_cov_from_ranks(_ranks, _g, 3, 8, "exact")
    st_matrix("V_dir", _V_dir)
}
matrix V_diff = V_cmd - V_dir
mata: st_numscalar("_r_vdiff", max(abs(st_matrix("V_diff"))))
test_assert "[identity] e(moment_vcov) equals direct vector jackknife covariance" scalar(_r_vdiff) < 1e-10

version 17
do "tests/test_framework.do"
test_section "PLQR finite score signs survive overflowed thresholds"
clear
set obs 6
gen double px=_n
gen double py=0
gen double pd=-1e100
gen double pzero=0
foreach sign in -1 1 {
    quietly replace pd=`sign'*1e100
    quietly plqr_dnn py, x(px) x0(0) d(pd) md(pzero) g(pzero) thetaeval(1e220) s(2)
    matrix PO_SCORE=e(model_signals)
    matrix PO_MOMENT=e(moment)
    mata: st_numscalar("PO_SIGN",all(st_matrix("PO_SCORE"):/(-5e99):==1))
    test_assert "`sign': one-dimensional overflow has the analytic represented score" PO_SIGN
    test_assert "`sign': overflow retains the finite localized moment" abs(PO_MOMENT[1,1]/(-5e99)-1)<1e-15
    quietly rdistnn_predict, x0(3) sefit
    matrix PO_FIT=r(fit)
    test_assert "`sign': frozen prediction retains the recovered score" abs(PO_FIT[1,1]/(-5e99)-1)<1e-15
}
quietly replace pd=-1e100
quietly plqr_dnn py, x(px) x0(0) d(pd) md(pzero) g(pzero) thetaeval(1e208) s(2)
matrix PO_SCORE=e(model_signals)
mata: st_numscalar("PO_SIGN",all(st_matrix("PO_SCORE"):/(-5e99):==1))
test_assert "finite negative threshold at -1e308 retains its existing indicator" PO_SIGN

gen double pd2=-1e100
quietly replace pd=1e100
quietly plqr_dnn py, x(px) x0(0) d(pd pd2) md(pzero pzero) g(pzero) thetaeval(.5e240 1e240) s(2)
matrix PO_SCORE=e(model_signals)
mata: st_numscalar("PO_SIGN",all(st_matrix("PO_SCORE")[,1]:/(5e99):==1) & all(st_matrix("PO_SCORE")[,2]:/(-5e99):==1))
test_assert "decisive mixed-sign overflow retains both analytic component signs" PO_SIGN

test_section "PLQR exact huge cancellation and input-error state"
gen double pg=-1
capture noisily plqr_dnn py, x(px) x0(0) d(pd pd2) md(pzero pzero) g(pg) thetaeval(1e240 1e240) s(2)
local exact_rc=_rc
test_assert "huge products retain the decisive unit residual" `exact_rc'==0
if `exact_rc'==0 {
    matrix PO_EXACT=e(model_signals)
    mata: st_numscalar("PO_EXACT_SIGN",all(st_matrix("PO_EXACT")[,1]:/(5e99):==1) & all(st_matrix("PO_EXACT")[,2]:/(-5e99):==1))
    test_assert "exact product cancellation produces the analytic finite scores" PO_EXACT_SIGN
}
do "tests/helpers/postestimation-state.do"
tempfile caller
quietly save "`caller'", replace
_pmpe_snapshot
capture noisily plqr_dnn py, x(px) x0(0) d(pd pd2) md(pzero pzero) g(pg) thetaeval(1) s(2)
local rc=_rc
test_assert "unaligned supplied comparison coefficients return an input error" `rc'==198
_pmpe_unchanged, data("`caller'")
local unchanged=r(ok)
test_assert "comparison input error preserves complete parent e() and caller data" `unchanged'
_rdistnn_test_e_snapshot_drop, name(PMPE_ESTATE)
macro drop PMPE_*

// The inclusive ordinary comparison is shared across all score owners.
test_section "PLQR ordinary equality"
clear
set obs 6
gen double pr_x=_n
gen double pr_d=1
gen double pr_zero=0
gen double pr_y=1
quietly plqr_dnn pr_y, x(pr_x) x0(0) d(pr_d) md(pr_zero) g(pr_zero) thetaeval(1) s(2)
matrix PR_SCORE=e(model_signals)
mata: st_numscalar("PR_ORDINARY",all(st_matrix("PR_SCORE"):==-.5))
test_assert "ordinary finite equality retains the existing inclusive indicator" PR_ORDINARY

test_section "PLQR exact predicate resolves ordinary cancellation"
clear
set obs 6
gen double ex_x=_n
gen double ex_y=.5
gen double ex_d1=18014398509481984
gen double ex_d2=1
gen double ex_d3=-ex_d1
gen double ex_zero=0
quietly plqr_dnn ex_y, x(ex_x) x0(0) d(ex_d1 ex_d2 ex_d3) md(ex_zero ex_zero ex_zero) g(ex_zero) thetaeval(1 1 1) s(2)
matrix EX_SCORE=e(model_signals)
// Exact dot product is1, margin1/2: the inclusive indicator is true.
mata: st_numscalar("EX_ORDINARY",all(st_matrix("EX_SCORE")[,1]:==-9007199254740992) & all(st_matrix("EX_SCORE")[,2]:==-.5) & all(st_matrix("EX_SCORE")[,3]:==9007199254740992))
test_assert "ordinary cancellation retains the exact half-unit comparison margin" EX_ORDINARY
quietly replace ex_d1=8e307
quietly replace ex_d2=-8e307
quietly replace ex_y=0
gen double ex_g=0
foreach residual in zero positive negative {
    quietly replace ex_g=0
    if "`residual'"=="positive" mata: st_store(.,"ex_g",J(st_nobs(),1,rdistnn_power2(-1074)))
    if "`residual'"=="negative" quietly replace ex_g=-1
    capture noisily plqr_dnn ex_y, x(ex_x) x0(0) d(ex_d1 ex_d2) md(ex_zero ex_zero) g(ex_g) thetaeval(1 1) s(2)
    local exact_rc=_rc
    test_assert "`residual': exact huge cancellation is supported" `exact_rc'==0
    if `exact_rc'==0 {
        local expected=-4e307
        if "`residual'"=="negative" local expected=4e307
        matrix EX_SCORE=e(model_signals)
        mata: st_numscalar("EX_SIGN",all(st_matrix("EX_SCORE")[,1]:/(`expected'):==1) & all(st_matrix("EX_SCORE")[,2]:/(-(`expected')):==1))
        test_assert "`residual': equality or decisive tiny residual gives exact finite scores" EX_SIGN
    }
}
test_section "PLQR exact predicate retains products below the double range"
foreach direction in -1 1 {
    capture noisily mata: EX_BELOW=rdistnn_plqr_indicator((8e307,8e307,rdistnn_power2(-1074)),(8e307 \ -8e307 \ (`direction'*rdistnn_power2(-1074))),0,0); st_numscalar("EX_BELOW",EX_BELOW)
    local exact_rc=_rc
    test_assert "`direction': product2^-2148 decides the exact comparison" `exact_rc'==0 & EX_BELOW==(`direction'>0)
}
