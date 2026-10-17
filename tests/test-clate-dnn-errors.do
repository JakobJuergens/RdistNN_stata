do "tests/test_framework.do"
clear
set obs 12
gen double x = _n
gen double Y = sin(_n)
gen byte A = mod(_n,2)
gen byte D = A
gen double e = .5
gen double y0 = 0
gen double y1 = 1
gen double d0 = 0
gen double d1 = 1
local base "xcols(x) x0(0) dcol(D) acol(A) ecol(e) my0col(y0) my1col(y1) md1col(d1) s(3)"
test_section "clate_dnn: validation"
capture clate_dnn Y, `base'
test_assert "two-sided needs zero-arm prediction" _rc==198
capture clate_dnn Y, `base' onesided md0col(d0) clip(.5)
test_assert "clip bounded away from half" _rc==198
capture clate_dnn Y, `base' onesided denominatortol(-1)
test_assert "negative tolerance rejected" _rc==198
capture clate_dnn Y, `base' onesided training(magic)
test_assert "unknown training label rejected" _rc==198
capture clate_dnn Y, `base' onesided pseudooutcome denominator(d1)
test_assert "mixed input routes rejected" _rc==198
capture clate_dnn Y, xcols(x) x0(0) pseudooutcome
test_assert "paired route needs denominator" _rc==198
capture clate_dnn Y, xcols(x) x0(0) denominator(d1)
test_assert "denominator requires paired route" _rc==198
replace D = 2 in 1
capture clate_dnn Y, `base' onesided
test_assert "nonbinary treatment rejected" _rc==198
replace D = 1 in 1
replace A = 0 in 1
capture clate_dnn Y, `base' onesided
test_assert "one-sided structural violation rejected" _rc==198
replace D = 0 in 1
replace e = 1.1 in 1
capture clate_dnn Y, `base' onesided
test_assert "unclipped propensity outside range rejected" _rc==198
replace e = . in 1
capture clate_dnn Y, `base' onesided
test_assert "missing propensity rejected" _rc==198
replace e = .5 in 1
replace d1 = 1.1 in 1
capture clate_dnn Y, `base' onesided
test_assert "treatment means outside range rejected" _rc==198
replace d1 = 1 in 1
replace d0 = .1 in 1
capture clate_dnn Y, `base' onesided md0col(d0)
test_assert "one-sided supplied zero-arm must be zero" _rc==198
replace d0 = 0 in 1
capture clate_dnn Y, `base' onesided x(x)
test_assert "covariate aliases cannot conflict" _rc==198
capture clate_dnn Y, `base' onesided x0matrix(MISSING)
test_assert "target specifications cannot conflict" _rc==198
capture clate_dnn Y, xcols(x) x0(0) pseudooutcome denominator(d1) s(12)
test_assert "inference requires s less than n" _rc==198
clate_dnn Y, xcols(x) x0(0) pseudooutcome denominator(d1) s(12) nose
test_assert "nose permits s equal n" e(s)==12
