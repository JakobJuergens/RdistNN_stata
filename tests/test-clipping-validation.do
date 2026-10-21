version 17
do "tests/helpers/fitted-data.do"
test_section "propensity clipping intervals must remain representable"
foreach family in cate clate {
    local extra "pseudooutcome"
    if "`family'"=="clate" local extra "`extra' denominator(fr_den)"
    capture quietly `family'_dnn fr_y, xcols(fr_x) x0(0) s(3) `extra' clip(1e-20)
    local rc=_rc
    test_assert "`family': unrepresentable clipping interval rejected" `rc'==198
}
