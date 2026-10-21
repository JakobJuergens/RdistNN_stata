version 17
do "tests/test_framework.do"
test_section "estimator and derived-fit stored-result metadata"
clear
set obs 20
gen double mx=(_n-10)/10
gen double my=1+mx^2+.1*sin(_n)
gen double md=.8+.05*mx
foreach family in np tdnn cate clate {
    local cmd "`family'_dnn"
    local extra "s(4)"
    local vce "rank_shift"
    local label "Jackknife"
    if "`family'"=="tdnn" {
        local cmd "np_tdnn"
        local extra "s1(2) s2(4)"
        local vce "signed_rank_shift"
    }
    if "`family'"=="cate" {
        local extra "`extra' pseudooutcome"
        local label "Frozen-score jackknife"
    }
    if "`family'"=="clate" {
        local extra "`extra' pseudooutcome denominator(md)"
        local vce "paired_rank_shift"
        local label "Paired frozen-score jackknife"
    }
    quietly `cmd' my, xcols(mx) x0(0) `extra'
    test_assert "`family': depvar records the user input" "`e(depvar)'"=="my"
    test_assert "`family': vce identifies the calculation" "`e(vce)'"=="`vce'"
    test_assert "`family': vcetype describes uncertainty" "`e(vcetype)'"=="`label'"
    quietly `cmd'
    test_assert "`family': replay preserves metadata" "`e(depvar)'"=="my" & "`e(vce)'"=="`vce'" & "`e(vcetype)'"=="`label'"
    quietly rdistnn_predict, x0(.2) sefit post
    test_assert "`family': sefit retains input and VCE" "`e(depvar)'"=="my" & "`e(vce)'"=="`vce'" & "`e(vcetype)'"=="`label'"
    quietly rdistnn_predict, x0(.1) post
    test_assert "`family': point-only derived fit records none" "`e(depvar)'"=="my" & "`e(vce)'"=="none" & "`e(vcetype)'"==""
    quietly rdistnn_predict, x0(.1) sefit post
    test_assert "`family': reactivating SE restores label" "`e(vce)'"=="`vce'" & "`e(vcetype)'"=="`label'"
    quietly `cmd' my, xcols(mx) x0(0) `extra' nose
    test_assert "`family': nose has outcome and no VCE label" "`e(depvar)'"=="my" & "`e(vce)'"=="none" & "`e(vcetype)'"==""
}
test_section "CLATE metadata with undefined ratios and augmented signals"
replace md=0
quietly clate_dnn my, denominator(md) pseudooutcome xcols(mx) x0(0) s(4)
test_assert "undefined ratio retains user input name" "`e(depvar)'"=="my"
test_assert "undefined ratio retains requested covariance method" "`e(vce)'"=="paired_rank_shift"
test_assert "method metadata does not override unavailable inference" missing(el(e(estimates),1,1)) & el(e(status),1,5)==0
gen byte ma=mod(_n,2)
gen double me=.5
gen double my0=mx
gen double my1=mx+1
gen double md0=0
gen double md1=1
quietly clate_dnn my, xcols(mx) x0(0) dcol(ma) acol(ma) ecol(me) ///
    my0col(my0) my1col(my1) md0col(md0) md1col(md1) s(4)
test_assert "augmented IV depvar identifies the outcome" "`e(depvar)'"=="my"
test_assert "augmented IV identifies paired frozen scores" "`e(vce)'"=="paired_rank_shift" & "`e(vcetype)'"=="Paired frozen-score jackknife"
