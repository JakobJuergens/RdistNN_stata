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
plqr_dnn y, x(x1 x2) x0(0 0) d(d1 d2) q(q) md(md1 md2) ///
    beta(1 -0.5) tau(0.5) k(4) s(12)

matrix _b = e(b)
matrix _V = e(V)
matrix _beta = e(beta)
local cmd "`e(cmd)'"
local kd = e(kd)
local tau = e(tau)
local nb = colsof(_b)
local bcols : colnames _b
local brows : rownames _V
local bcolsV : colnames _V

scalar _b1 = _b[1,1]
scalar _b2 = _b[1,2]
scalar _v11 = _V[1,1]
scalar _v22 = _V[2,2]
scalar _v12 = _V[1,2]
scalar _v21 = _V[2,1]
scalar _beta1 = _beta[1,1]
scalar _beta2 = _beta[1,2]

test_assert "[smoke] e(cmd) = plqr_dnn" "`cmd'" == "plqr_dnn"
test_assert "[smoke] e(b) has two moment components" `nb' == 2
test_assert "[smoke] e(b) column names follow d()" "`bcols'" == "d1 d2"
test_assert "[smoke] e(V) row names follow d()" "`brows'" == "d1 d2"
test_assert "[smoke] e(V) column names follow d()" "`bcolsV'" == "d1 d2"
test_assert "[smoke] moments are finite" scalar(_b1) < . & scalar(_b2) < .
test_assert "[smoke] variances are finite and non-negative" ///
    scalar(_v11) < . & scalar(_v11) >= 0 & scalar(_v22) < . & scalar(_v22) >= 0
test_assert "[smoke] covariance matrix is symmetric" abs(scalar(_v12) - scalar(_v21)) < 1e-10
test_assert "[smoke] e(beta) stores supplied beta" ///
    abs(scalar(_beta1) - 1) < 1e-12 & abs(scalar(_beta2) + 0.5) < 1e-12
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

plqr_dnn y, x(x1 x2) x0(0 0) d(d1 d2) q(q) md(md1 md2) ///
    beta(1 -0.5) tau(0.5) k(4) s(10)
matrix _b = e(b)
scalar _r1 = abs(_b[1,1] - `m1')
scalar _r2 = abs(_b[1,2] - `m2')
test_assert "[identity] first moment equals manual weighted sum" scalar(_r1) < 1e-8
test_assert "[identity] second moment equals manual weighted sum" scalar(_r2) < 1e-8

// ════════════════════════════════════════════════════════════
test_section "plqr_dnn: [identity] covariance equals direct Mata helper"
// ════════════════════════════════════════════════════════════

rdnn_dgp_plqr, n(80) seed(403)
plqr_dnn y, x(x1 x2) x0(0.2 -0.1) d(d1 d2) q(q) md(md1 md2) ///
    beta(1 -0.5) tau(0.5) k(3) s(8)
matrix _V_cmd = e(V)

mata {
    _X     = st_data(., "x1 x2")
    _D     = st_data(., "d1 d2")
    _md    = st_data(., "md1 md2")
    _q     = st_data(., "q")
    _y     = st_data(., "y")
    _beta  = (1 \ -0.5)
    _x0    = (0.2, -0.1)
    _ranks = rdistnn_ranks(_X, _x0, "l2")
    _resid = 0.5 :- (_y :<= (_D * _beta :+ _q))
    _g     = (_D :- _md) :* (_resid * J(1, cols(_D), 1))
    _V_dir = rdistnn_jk_cov_from_ranks(_ranks, _g, 3, 8, "exact")
    st_matrix("_V_dir", _V_dir)
}
matrix _V_diff = _V_cmd - _V_dir
mata: st_numscalar("_r_vdiff", max(abs(st_matrix("_V_diff"))))
test_assert "[identity] e(V) equals direct vector jackknife covariance" scalar(_r_vdiff) < 1e-10
