// test-plqr-dnn-r-identities.do
// Additional PLQR/Z-moment checks ported from the R score and user tests.

quietly do "src/_rdistnn_mata.ado"
do "tests/test_framework.do"

test_section "plqr_dnn: [edge] p=1, kd=1 command shape and identity"

clear
set obs 60
set seed 1101
gen double x1 = rnormal()
gen double md1 = 0.2 * x1
gen double d1 = md1 + rnormal() * 0.3
gen double q = 0.1 * x1
gen double y = 0.75 * d1 + q + rnormal() * 0.4

weights_dnn x1, x0(-0.1) k(2) s(6) generate(w_plqr_1d)
gen double g1 = (0.5 - (y <= 0.75 * d1 + q)) * (d1 - md1)
gen double wg1 = w_plqr_1d * g1
quietly summarize wg1, meanonly
local m1 = r(sum)

plqr_dnn y, x(x1) x0(-0.1) d(d1) g(q) md(md1) thetaeval(0.75) k(2) s(6)
matrix B1 = e(moment)
matrix V1 = e(moment_vcov)
local kd1 = e(kd)
local kx1 = e(kx)
test_assert "[edge] one-dimensional PLQR moment equals manual sum" abs(B1[1,1] - `m1') < 1e-8
test_assert "[edge] one-dimensional PLQR e(kx)=1" `kx1' == 1
test_assert "[edge] one-dimensional PLQR e(kd)=1" `kd1' == 1
test_assert "[edge] one-dimensional PLQR variance finite" V1[1,1] < . & V1[1,1] >= 0

test_section "plqr_dnn: [indicator] threshold equality uses <= branch"

clear
set obs 30
gen double x1 = _n / 30
gen double x2 = 1 - x1
gen double d1 = 1 + x1
gen double d2 = -0.5 + x2
gen double md1 = 0
gen double md2 = 0
gen double q = 0.25
gen double y = d1 - 0.5 * d2 + q

weights_dnn x1 x2, x0(0.4 0.6) k(2) s(5) generate(w_eq)
gen double g_eq1 = -0.5 * d1
gen double g_eq2 = -0.5 * d2
gen double wg_eq1 = w_eq * g_eq1
gen double wg_eq2 = w_eq * g_eq2
quietly summarize wg_eq1, meanonly
local m_eq1 = r(sum)
quietly summarize wg_eq2, meanonly
local m_eq2 = r(sum)

plqr_dnn y, x(x1 x2) x0(0.4 0.6) d(d1 d2) g(q) md(md1 md2) ///
    thetaeval(1 -0.5) tau(0.5) k(2) s(5)
matrix Beq = e(moment)
test_assert "[indicator] equality contributes tau-1 for d1" abs(Beq[1,1] - `m_eq1') < 1e-8
test_assert "[indicator] equality contributes tau-1 for d2" abs(Beq[1,2] - `m_eq2') < 1e-8

test_section "plqr_dnn: [edge] zero residualized regressors give zero moment and V"

rdnn_dgp_plqr, n(75) seed(1102)
replace d1 = md1
replace d2 = md2
plqr_dnn y, x(x1 x2) x0(0 0) d(d1 d2) g(q) md(md1 md2) ///
    thetaeval(1 -0.5) tau(0.5) k(2) s(7)
matrix Bzero = e(moment)
matrix Vzero = e(moment_vcov)
mata {
    st_numscalar("_plqr_zero_b", max(abs(st_matrix("Bzero"))))
    st_numscalar("_plqr_zero_v", max(abs(st_matrix("Vzero"))))
}
test_assert "[edge] zero D-md gives zero PLQR moments" scalar(_plqr_zero_b) < 1e-12
test_assert "[edge] zero D-md gives zero PLQR covariance" scalar(_plqr_zero_v) < 1e-12

test_section "plqr_dnn: [theta] changing theta_eval can change moments"

rdnn_dgp_plqr, n(100) seed(1103)
plqr_dnn y, x(x1 x2) x0(0.1 -0.2) d(d1 d2) g(q) md(md1 md2) ///
    thetaeval(1 -0.5) tau(0.5) k(3) s(9)
matrix Mtheta1 = e(moment)
plqr_dnn y, x(x1 x2) x0(0.1 -0.2) d(d1 d2) g(q) md(md1 md2) ///
    thetaeval(0.2 0.2) tau(0.5) k(3) s(9)
matrix Mtheta2 = e(moment)
matrix MthetaDiff = Mtheta1 - Mtheta2
mata: st_numscalar("_plqr_theta_diff", max(abs(st_matrix("MthetaDiff"))))
test_assert "[theta] supplied theta_eval enters the PLQR indicator" scalar(_plqr_theta_diff) > 1e-8

test_section "plqr_dnn: [markout] missing y/D/g/m_D rows are excluded"

rdnn_dgp_plqr, n(90) seed(1104)
replace y = . in 2
replace d1 = . in 3
replace q = . in 4
replace md2 = . in 5
plqr_dnn y, x(x1 x2) x0(0 0) d(d1 d2) g(q) md(md1 md2) ///
    thetaeval(1 -0.5) tau(0.5) k(2) s(6)
local n_mark = e(N)
quietly count if e(sample)
local n_esample = r(N)
test_assert "[markout] missing PLQR source rows are dropped" `n_mark' == 86
test_assert "[markout] PLQR e(sample) agrees with e(N)" `n_esample' == 86
