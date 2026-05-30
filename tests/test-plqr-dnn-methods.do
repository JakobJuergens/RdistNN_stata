// test-plqr-dnn-methods.do
// Tests for plqr_dnn method, distance, and sample-restriction options.
//
// Tags: [plqr] [opts] [distance] [ifin]

quietly do "src/_rdistnn_mata.ado"
do "tests/test_framework.do"

// ════════════════════════════════════════════════════════════
test_section "plqr_dnn: [opts] exact and asymp methods"
// ════════════════════════════════════════════════════════════

rdnn_dgp_plqr, n(130) seed(501)
plqr_dnn y, x(x1 x2) x0(0 0) d(d1 d2) q(q) md(md1 md2) ///
    beta(1 -0.5) tau(0.5) k(4) s(12) method(exact)
matrix _b_exact = e(b)
local meth_exact "`e(method)'"

plqr_dnn y, x(x1 x2) x0(0 0) d(d1 d2) q(q) md(md1 md2) ///
    beta(1 -0.5) tau(0.5) k(4) s(12) method(asymp)
matrix _b_asymp = e(b)
local meth_asymp "`e(method)'"

scalar _exact_fin = (_b_exact[1,1] < . & _b_exact[1,2] < .)
scalar _asymp_fin = (_b_asymp[1,1] < . & _b_asymp[1,2] < .)

test_assert "[opts] exact method recorded" "`meth_exact'" == "exact"
test_assert "[opts] asymp method recorded" "`meth_asymp'" == "asymp"
test_assert "[opts] exact moments finite" scalar(_exact_fin) == 1
test_assert "[opts] asymp moments finite" scalar(_asymp_fin) == 1

// ════════════════════════════════════════════════════════════
test_section "plqr_dnn: [distance] l2 and mahalanobis"
// ════════════════════════════════════════════════════════════

rdnn_dgp_plqr, n(130) seed(502)
plqr_dnn y, x(x1 x2) x0(0.1 -0.2) d(d1 d2) q(q) md(md1 md2) ///
    beta(1 -0.5) tau(0.5) k(4) s(12) distance(l2)
matrix _b_l2 = e(b)
local dist_l2 "`e(distance)'"

plqr_dnn y, x(x1 x2) x0(0.1 -0.2) d(d1 d2) q(q) md(md1 md2) ///
    beta(1 -0.5) tau(0.5) k(4) s(12) distance(mahalanobis)
matrix _b_mh = e(b)
local dist_mh "`e(distance)'"

scalar _l2_fin = (_b_l2[1,1] < . & _b_l2[1,2] < .)
scalar _mh_fin = (_b_mh[1,1] < . & _b_mh[1,2] < .)

test_assert "[distance] l2 recorded" "`dist_l2'" == "l2"
test_assert "[distance] mahalanobis recorded" "`dist_mh'" == "mahalanobis"
test_assert "[distance] l2 moments finite" scalar(_l2_fin) == 1
test_assert "[distance] mahalanobis moments finite" scalar(_mh_fin) == 1

// ════════════════════════════════════════════════════════════
test_section "plqr_dnn: [ifin] restricted sample matches manual weights"
// ════════════════════════════════════════════════════════════

rdnn_dgp_plqr, n(150) seed(503)
gen byte keep_plqr = (_n <= 120 & x1 > -0.25)

weights_dnn x1 x2 if keep_plqr, x0(0.3 -0.1) k(3) s(9) generate(w_plqr_restricted)
gen double g1 = (0.5 - (y <= d1 - 0.5 * d2 + q)) * (d1 - md1) if keep_plqr
gen double g2 = (0.5 - (y <= d1 - 0.5 * d2 + q)) * (d2 - md2) if keep_plqr
gen double wg1 = w_plqr_restricted * g1 if keep_plqr
gen double wg2 = w_plqr_restricted * g2 if keep_plqr
quietly summarize wg1 if keep_plqr, meanonly
local m1 = r(sum)
quietly summarize wg2 if keep_plqr, meanonly
local m2 = r(sum)

plqr_dnn y if keep_plqr, x(x1 x2) x0(0.3 -0.1) d(d1 d2) q(q) md(md1 md2) ///
    beta(1 -0.5) tau(0.5) k(3) s(9)
matrix _b = e(b)
scalar _r1 = abs(_b[1,1] - `m1')
scalar _r2 = abs(_b[1,2] - `m2')

test_assert "[ifin] first restricted moment matches manual sum" scalar(_r1) < 1e-8
test_assert "[ifin] second restricted moment matches manual sum" scalar(_r2) < 1e-8
