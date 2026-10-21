// test-np-dnn-r-identities.do
// Additional np_dnn command identities inspired by the R user-facing tests.

quietly do "src/_rdistnn_mata.ado"
do "tests/test_framework.do"

test_section "np_dnn: [identity] weighted-mean identity across methods/distances"

rdnn_dgp_np, n(140) seed(901)

weights_dnn x1 x2, x0(0.25 -0.4) k(3) s(9) method(asymp) generate(w_np_ap)
mata {
    _w = st_data(., "w_np_ap")
    _y = st_data(., "y")
    st_numscalar("_np_manual_ap", sum(_w :* _y))
}
np_dnn y, xcols(x1 x2) x0(0.25 -0.4) k(3) s(9) method(asymp)
local mu_ap = _b[mu]
test_assert "[identity] asymp mu equals sum(w*y)" abs(`mu_ap' - scalar(_np_manual_ap)) < 1e-8

weights_dnn x1 x2, x0(0.25 -0.4) k(3) s(9) distance(mahalanobis) generate(w_np_mh)
mata {
    _w = st_data(., "w_np_mh")
    _y = st_data(., "y")
    st_numscalar("_np_manual_mh", sum(_w :* _y))
}
np_dnn y, xcols(x1 x2) x0(0.25 -0.4) k(3) s(9) distance(mahalanobis)
local mu_mh = _b[mu]
test_assert "[identity] mahalanobis mu equals sum(w*y)" abs(`mu_mh' - scalar(_np_manual_mh)) < 1e-8

test_section "np_dnn: [ifin] restricted sample identity and e(sample)"

rdnn_dgp_np, n(130) seed(902)
gen byte keep_np = (_n >= 10 & _n <= 100 & x2 < 0.75)
quietly count if keep_np
local n_keep = r(N)

weights_dnn x1 x2 if keep_np, x0(-0.1 0.2) k(2) s(7) generate(w_np_if)
mata {
    _w = st_data(., "w_np_if", "keep_np")
    _y = st_data(., "y", "keep_np")
    st_numscalar("_np_manual_if", sum(_w :* _y))
}
np_dnn y if keep_np, xcols(x1 x2) x0(-0.1 0.2) k(2) s(7)
local mu_if = _b[mu]
local n_if = e(N)
quietly count if e(sample)
local n_esample = r(N)
test_assert "[ifin] restricted mu equals manual restricted weighted mean" abs(`mu_if' - scalar(_np_manual_if)) < 1e-8
test_assert "[ifin] e(N) matches restricted sample" `n_if' == `n_keep'
test_assert "[ifin] e(sample) marks exactly restricted sample" `n_esample' == `n_keep'

test_section "np_dnn: [markout] missing outcomes/covariates are excluded"

rdnn_dgp_np, n(80) seed(903)
replace y = . in 3
replace x1 = . in 4
replace x2 = . in 5
np_dnn y, xcols(x1 x2) x0(0 0) k(2) s(6)
local n_mark = e(N)
quietly count if e(sample)
local n_sample_mark = r(N)
test_assert "[markout] missing y/x rows are dropped from e(N)" `n_mark' == 77
test_assert "[markout] e(sample) agrees with e(N) after markout" `n_sample_mark' == 77

test_section "np_dnn: [edge] singleton k=s=1 uniform U-statistic identity"

clear
set obs 6
gen double x1 = _n
gen double y = 10 + _n
np_dnn y, xcols(x1) x0(3.2) k(1) s(1)
local nn_mu = _b[mu]
local nn_se = _se[mu]
test_assert "[edge] k=s=1 returns the sample mean" abs(`nn_mu' - 13.5) < 1e-12
test_assert "[edge] k=s=1 standard error is finite" `nn_se' < .
