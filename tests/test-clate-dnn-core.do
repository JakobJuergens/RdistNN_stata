// Independent four-row ratio and perfect-compliance reductions.
do "tests/test_framework.do"
clear
set obs 4
gen double x = _n
gen double py = cond(_n==1,2,cond(_n==2,1,cond(_n==3,4,3)))
gen double pd = cond(_n==2,-1,cond(_n==3,2,1))

test_section "clate_dnn: paired analytic fixture"
clate_dnn py, denominator(pd) pseudooutcome xcols(x) x0(0) s(2)
test_assert "ratio is four" abs(_b[theta]-4)<1e-12
test_assert "paired variance is 34" abs(el(e(V),1,1)-34)<1e-10
test_assert "reduced form is two" abs(el(e(components),1,1)-2)<1e-12
test_assert "first stage is one half" abs(el(e(components),1,2)-.5)<1e-12
test_assert "own estimator identity" "`e(cmd)'"=="clate_dnn"
test_assert "paired covariance declaration" "`e(se_method)'"=="paired_rank_shift"
test_assert "all observations in sample" e(N)==4
clate_dnn, level(90)

matrix CT = (0 \ 5)
clate_dnn py, denominator(pd) pseudooutcome xcols(x) x0matrix(CT) s(2)
matrix CB = e(estimates)
matrix CV = e(theta_vcov)
test_assert "two target coefficients" colsof(e(b))==2
test_assert "batch covariance symmetric" abs(el(CV,1,2)-el(CV,2,1))<1e-12
clate_dnn py, denominator(pd) pseudooutcome xcols(x) x0(5) s(2)
test_assert "batch second estimate equals scalar" abs(_b[theta]-el(CB,1,2))<1e-12
test_assert "batch second variance equals scalar" abs(el(e(V),1,1)-el(CV,2,2))<1e-12

test_section "clate_dnn: signed and zero denominators"
replace pd = -pd
clate_dnn py, denominator(pd) pseudooutcome xcols(x) x0(0) s(2)
test_assert "negative stage retained" abs(_b[theta]+4)<1e-12
test_assert "negative stage flagged" el(e(status),1,2)==1
replace pd = 0
clate_dnn py, denominator(pd) pseudooutcome xcols(x) x0(0) s(2)
test_assert "zero stage ratio unavailable" missing(el(e(estimates),1,1))
test_assert "zero stage status" el(e(status),1,1)==1
test_assert "zero stage components retained" abs(el(e(components),1,1)-2)<1e-12
test_assert "zero stage has no Wald inference" el(e(status),1,5)==0
clate_dnn

test_section "clate_dnn: perfect compliance reduces to CATE"
clear
set obs 40
gen double x = (_n-20)/10
gen byte A = mod(_n,2)
gen byte D = A
gen double e = .4
gen double y0 = x/3
gen double y1 = y0 + 1 + sin(x)
gen double Y = y0 + D*(y1-y0) + cos(_n)/5
gen double d0 = 0
gen double d1 = 1
cate_dnn Y, xcols(x) x0(0) tcol(A) ecol(e) m0col(y0) m1col(y1) s(5)
scalar cate_b = _b[tau]
scalar cate_v = el(e(V),1,1)
clate_dnn Y, xcols(x) x0(0) dcol(D) acol(A) ecol(e) ///
    my0col(y0) my1col(y1) md0col(d0) md1col(d1) s(5)
test_assert "perfect compliance CATE estimate" abs(_b[theta]-scalar(cate_b))<1e-12
test_assert "perfect compliance CATE variance" abs(el(e(V),1,1)-scalar(cate_v))<1e-12
test_assert "unit first stage" abs(el(e(components),1,2)-1)<1e-12
test_assert "instrument overlap available" el(e(diagnostics),1,12)==1
test_assert "instrument label" "`e(overlap_role)'"=="instrument"
matrix CLB = e(b)
matrix CLV = e(V)
estimates store CLATE_saved
cate_dnn Y, xcols(x) x0(0) tcol(A) ecol(e) m0col(y0) m1col(y1) s(5)
estimates restore CLATE_saved
clate_dnn
test_assert "restored coefficient unchanged" abs(_b[theta]-el(CLB,1,1))<1e-12

clate_dnn Y if _n<=30, xcols(x) x0(0) dcol(D) acol(A) ecol(e) ///
    my0col(y0) my1col(y1) md1col(d1) onesided s(5) nose
test_assert "if restriction retained" e(N)==30
test_assert "onesided metadata" e(one_sided)==1
test_assert "nose metadata" "`e(jackknife)'"=="none"
test_assert "nose standard error missing" missing(el(e(se_targets),1,1))

test_section "clate_dnn: partially undefined batch keeps target identities"
clear
set obs 4
gen double x = _n
gen double py = 1
gen double pd = cond(_n<=2,0,1)
matrix CT = (0 \ 5)
clate_dnn py, denominator(pd) pseudooutcome xcols(x) x0matrix(CT) s(4) nose
test_assert "first undefined target retained" missing(el(e(estimates),1,1))
test_assert "second target remains defined" abs(el(e(estimates),1,2)-1)<1e-12
test_assert "only defined target posted" colsof(e(b))==1
test_assert "posted coefficient has target two identity" abs(_b[theta2]-1)<1e-12

replace py = cond(_n==1,2,cond(_n==2,1,cond(_n==3,4,3)))
replace pd = cond(_n<=2,0,cond(_n==3,-1,2))
clate_dnn py, denominator(pd) pseudooutcome xcols(x) x0matrix(CT) s(3)
test_assert "partial covariance fit posts defined target only" colsof(e(b))==1 & colsof(e(V))==1
test_assert "partial covariance fit preserves target two" abs(_b[theta2]-2.6)<1e-12
test_assert "partial covariance uses canonical target two variance" abs(el(e(V),1,1)-el(e(theta_vcov),2,2))<1e-12
test_assert "undefined covariance row remains missing" missing(el(e(theta_vcov),1,1)) & missing(el(e(theta_vcov),1,2))

test_section "clate_dnn: finite ratio survives unavailable component covariance"
replace py = cond(_n==1,2,cond(_n==2,1,cond(_n==3,4,3)))*1e160
replace pd = 1
capture noisily clate_dnn py, denominator(pd) pseudooutcome xcols(x) x0(0) s(2)
local extreme_rc = _rc
test_assert "finite extreme ratio remains a successful fit" `extreme_rc'==0
if `extreme_rc'==0 {
    test_assert "extreme ratio matches finite R reference" abs(_b[theta]/1e160-2)<1e-12
    test_assert "extreme component covariance flagged unavailable" el(e(status),1,4)==3
    test_assert "extreme fit retains unavailable SE" missing(el(e(se_targets),1,1))
}
