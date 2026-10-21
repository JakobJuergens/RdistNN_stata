// test-cate-dnn-score-inputs.do
// Parity tests for the R cate_dnn score interfaces and estimate-only path.

quietly do "src/_rdistnn_mata.ado"
do "tests/test_framework.do"

test_section "cate_dnn: [score] precomputed AIPW pseudo-outcome"

rdnn_dgp_basic, n(120) seed(2101)
gen double eclip_score = min(max(e, 1e-8), 1 - 1e-8)
gen double rho_aipw = D * (y - m1) / eclip_score ///
    - (1-D) * (y - m0) / (1 - eclip_score) + (m1 - m0)

cate_dnn y, tcol(D) xcols(x1 x2) x0(0.2 -0.1) k(2) s(8) ///
    ecol(e) m0col(m0) m1col(m1)
local tau_aipw = _b[tau]
local se_aipw = _se[tau]
local score_aipw "`e(score_method)'"
test_assert "[meta] AIPW depvar identifies the outcome" "`e(depvar)'" == "y"
test_assert "[meta] AIPW VCE identifies rank-shift deletions" "`e(vce)'" == "rank_shift"
test_assert "[meta] AIPW VCE labels frozen scores" "`e(vcetype)'" == "Frozen-score jackknife"

cate_dnn rho_aipw, xcols(x1 x2) x0(0.2 -0.1) k(2) s(8) pseudooutcome
local tau_pseudo = _b[tau]
local se_pseudo = _se[tau]
local score_pseudo "`e(score_method)'"
test_assert "[meta] pseudo-outcome depvar identifies the supplied signal" "`e(depvar)'" == "rho_aipw"
test_assert "[meta] pseudo-outcome VCE labels frozen scores" "`e(vcetype)'" == "Frozen-score jackknife"

test_assert "[score] AIPW and pseudo-outcome estimates agree" abs(`tau_aipw' - `tau_pseudo') < 1e-10
test_assert "[score] AIPW and pseudo-outcome SEs agree" abs(`se_aipw' - `se_pseudo') < 1e-10
test_assert "[score] AIPW metadata recorded" "`score_aipw'" == "aipw"
test_assert "[score] pseudo-outcome metadata recorded" "`score_pseudo'" == "pseudo_outcome"

test_section "cate_dnn: [nose] estimate-only fit and R default k"

rdnn_dgp_basic, n(50) seed(2103)
gen double rho_nose = 2 + 0.5 * x1
gen double d2_nose = x1^2 + x2^2
sort d2_nose
local rho_nearest = rho_nose[1]

cate_dnn rho_nose, xcols(x1 x2) x0(0 0) s(50) pseudooutcome nose
local tau_nose = _b[tau]
local k_default = e(k)
local s_nose = e(s)
local se_method "`e(se_method)'"
capture matrix list e(V)
local rc_v = _rc

test_assert "[nose] s=n and k=1 selects the nearest pseudo-outcome" abs(`tau_nose' - `rho_nearest') < 1e-10
test_assert "[nose] default k matches R paper default" `k_default' == 1
test_assert "[nose] s=n is retained" `s_nose' == 50
test_assert "[nose] e(se_method) records none" "`se_method'" == "none"
test_assert "[nose] e(V) is not posted" `rc_v' != 0
test_assert "[meta] nose retains the supplied signal name" "`e(depvar)'" == "rho_nose"
test_assert "[meta] nose records no VCE" "`e(vce)'" == "none"
test_assert "[meta] nose has no uncertainty label" "`e(vcetype)'" == ""

test_section "cate_dnn: [markout] pseudo-outcome sample"

rdnn_dgp_basic, n(40) seed(2104)
gen double rho_mark = y
replace rho_mark = . in 2/4
replace x1 = . in 5
cate_dnn rho_mark, xcols(x1 x2) x0(0 0) k(1) s(5) pseudooutcome nose
local n_mark = e(N)
quietly count if e(sample)
local n_esample = r(N)
test_assert "[markout] missing X/pseudo-outcome rows are dropped" `n_mark' == 36
test_assert "[markout] e(sample) matches pseudo-outcome sample" `n_esample' == 36
