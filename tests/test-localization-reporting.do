version 17
do "tests/test_framework.do"
test_section "Localization reporting and nonuniform tied-distance reference"
preserve
clear
set obs 4
gen double x = (_n > 2)
gen double signal = cond(_n == 1, 1, cond(_n == 2, 3, cond(_n == 3, 10, 20)))
foreach command in np_dnn cate_dnn {
    local scoreopt ""
    if "`command'" == "cate_dnn" local scoreopt "pseudooutcome"
    quietly `command' signal, xcols(x) x0(0) k(1) s(2) `scoreopt'
    tempname original covariance
    matrix `original' = e(b)
    matrix `covariance' = e(V)
    // Full block weights are (5/12,5/12,1/12,1/12); literal deletions
    // are (7,17/3,2,2), giving variance (3/4)*sum((L-mean(L))^2)=59/4.
    test_assert "[ties] `command' exact block-averaged mean" abs(el(`original',1,1)-25/6)<1e-12
    test_assert "[ties] `command' explicit deletion covariance" abs(el(`covariance',1,1)-59/4)<1e-12
    tempfile report
    quietly log using "`report'", text replace name(stage1_report)
    _rdistnn_report_header
    quietly log close stage1_report
    tempname handle
    file open `handle' using "`report'", read text
    file read `handle' line
    local found 0
    while r(eof) == 0 {
        if strpos(`"`line'"', "s/n") local found 1
        file read `handle' line
    }
    file close `handle'
    test_assert "[report] `command' shows s/n" `found' == 1
    quietly replace signal = 3 in 1
    quietly replace signal = 1 in 2
    quietly `command' signal, xcols(x) x0(0) k(1) s(2) `scoreopt'
    test_assert "[ties] `command' identical-X permutation invariance" abs(el(e(b),1,1)-el(`original',1,1))<1e-12
    quietly replace signal = 1 in 1
    quietly replace signal = 3 in 2
}
restore
