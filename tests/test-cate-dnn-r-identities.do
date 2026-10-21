// test-cate-dnn-r-identities.do
// Additional cate_dnn identities inspired by the R CATE user/core tests.

quietly do "src/_rdistnn_mata.ado"
do "tests/test_framework.do"

test_section "cate_dnn: [identity] AIPW clipping identity at propensity boundaries"

rdnn_dgp_basic, n(120) seed(1001)
replace e = 0 in 1/3
replace e = 1 in 4/6
gen double eclip_r = min(max(e, 1e-8), 1 - 1e-8)
gen double phi_r = D * (y - m1) / eclip_r ///
    - (1-D) * (y - m0) / (1 - eclip_r) ///
    + (m1 - m0)

weights_dnn x1 x2, x0(0.2 -0.1) k(3) s(9) generate(w_cate_r)
gen double wphi_r = w_cate_r * phi_r
quietly summarize wphi_r, meanonly
local tau_manual = r(sum)

cate_dnn y, tcol(D) xcols(x1 x2) x0(0.2 -0.1) k(3) s(9) ///
    ecol(e) m0col(m0) m1col(m1)
local tau_cmd = _b[tau]
local se_cmd = _se[tau]
test_assert "[identity] tau uses clipped AIPW pseudo-outcome" ///
    abs(`tau_cmd' - `tau_manual') / max(1, abs(`tau_manual')) < 1e-12
test_assert "[identity] boundary-propensity SE is finite" `se_cmd' < .

test_section "cate_dnn: [ifin] restricted sample identity"

rdnn_dgp_basic, n(140) seed(1002)
gen byte keep_cate = (_n <= 110 & x1 > -0.5)
quietly count if keep_cate
local n_keep = r(N)

gen double eclip_c = min(max(e, 1e-8), 1 - 1e-8) if keep_cate
gen double phi_c = D * (y - m1) / eclip_c ///
    - (1-D) * (y - m0) / (1 - eclip_c) ///
    + (m1 - m0) if keep_cate
weights_dnn x1 x2 if keep_cate, x0(-0.2 0.15) k(2) s(7) generate(w_cate_if)
gen double wphi_c = w_cate_if * phi_c if keep_cate
quietly summarize wphi_c if keep_cate, meanonly
local tau_if_manual = r(sum)

cate_dnn y if keep_cate, tcol(D) xcols(x1 x2) x0(-0.2 0.15) k(2) s(7) ///
    ecol(e) m0col(m0) m1col(m1)
local tau_if_cmd = _b[tau]
local n_if = e(N)
quietly count if e(sample)
local n_esample = r(N)
test_assert "[ifin] restricted CATE equals manual restricted AIPW mean" abs(`tau_if_cmd' - `tau_if_manual') < 1e-8
test_assert "[ifin] e(N) matches restricted CATE sample" `n_if' == `n_keep'
test_assert "[ifin] e(sample) marks restricted CATE sample" `n_esample' == `n_keep'

test_section "cate_dnn: [ties] deterministic tied covariates give fixed tau"

clear
set obs 12
gen double x1 = 0
gen double x2 = 0
gen byte D = mod(_n, 2)
gen double y = _n
gen double e = 0.5
gen double m0 = y - 1
gen double m1 = y + 1

gen double eclip_t = min(max(e, 1e-8), 1 - 1e-8)
gen double phi_t = D * (y - m1) / eclip_t ///
    - (1-D) * (y - m0) / (1 - eclip_t) ///
    + (m1 - m0)
weights_dnn x1 x2, x0(0 0) k(2) s(4) generate(w_tie_cate)
gen double wphi_t = w_tie_cate * phi_t
quietly summarize wphi_t, meanonly
local tau_tie_manual = r(sum)

cate_dnn y, tcol(D) xcols(x1 x2) x0(0 0) k(2) s(4) ecol(e) m0col(m0) m1col(m1)
local tau_tie_cmd = _b[tau]
test_assert "[ties] tied-X CATE uses deterministic row-index weights" abs(`tau_tie_cmd' - `tau_tie_manual') < 1e-8

test_section "cate_dnn: [markout] missing propensity errors; other nuisance rows are excluded"

rdnn_dgp_basic, n(90) seed(1003)
replace y = . in 2
replace D = . in 3
replace e = . in 4
replace m0 = . in 5
replace m1 = . in 6
capture noisily cate_dnn y, tcol(D) xcols(x1 x2) x0(0 0) k(2) s(6) ecol(e) m0col(m0) m1col(m1)
test_assert "[markout] missing raw propensity is rejected" _rc == 198
replace e = .5 in 4
cate_dnn y, tcol(D) xcols(x1 x2) x0(0 0) k(2) s(6) ecol(e) m0col(m0) m1col(m1)
local n_mark = e(N)
quietly count if e(sample)
local n_esample_mark = r(N)
test_assert "[markout] missing y/D/m0/m1 rows are dropped" `n_mark' == 86
test_assert "[markout] e(sample) agrees with CATE markout count" `n_esample_mark' == 86
