// Supported nuisance construction and removed-option regression checks.
do "tests/test_framework.do"
test_section "cate_dnn: supported CATE score construction"
rdnn_dgp_basic, n(40) seed(2321)
foreach opt in mcol mhat {
    capture cate_dnn y, xcols(x1 x2) x0(0 0) tcol(D) ecol(e) `opt'(m0) s(5)
    test_assert "[score] pooled-mean option `opt' is rejected" _rc == 198
    capture cate_dnn y, xcols(x1 x2) x0(0 0) pseudooutcome `opt'(m0) s(5)
    test_assert "[score] supplied signals cannot bypass rejection of `opt'" _rc == 198
}
quietly cate_dnn y, xcols(x1 x2) x0(0 0) tcol(D) ecol(e) m0col(m0) m1col(m1) s(5)
test_assert "[score] supported nuisance route records AIPW" "`e(score_method)'" == "aipw"
test_assert "[score] pooled nuisance metadata is absent" "`e(m_col)'" == ""
capture program drop _cate_test_retired
program define _cate_test_retired, eclass
    ereturn scalar model_version = 0
    ereturn local score_method "retired_score"
end
_cate_test_retired
program drop _cate_test_retired
capture cate_dnn_plot, type(influence) nodraw
test_assert "[score] retired saved score cannot reconstruct influence" _rc == 498

preserve
clear
set obs 4
gen double x = _n
gen byte D = _n == 4
gen double y = D
gen double m0 = 0
gen double m1 = 1
gen double e = .24
foreach propensity in .24 .26 .30 {
    quietly replace e = `propensity'
    quietly cate_dnn y, xcols(x) x0(0) tcol(D) ecol(e) m0col(m0) m1col(m1) s(1) nose
    test_assert "[score] AIPW population identity with propensity `propensity'" abs(_b[tau]-1) < 1e-15
}
restore
