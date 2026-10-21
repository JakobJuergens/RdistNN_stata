// Independent concentration references for finite exceptional score ranges.
version 17
do "tests/test_framework.do"
do "tests/helpers/postestimation-state.do"
test_section "finite diagnostic contributions use one current policy"
foreach family in cate clate {
    clear
    set obs 12
    gen double range_x=_n
    gen double range_signal=cond(_n==11,8e307,-6e307)
    gen double range_den=1
    local extra "pseudooutcome"
    if "`family'"=="clate" local extra "`extra' denominator(range_den)"
    quietly `family'_dnn range_signal, `extra' xcols(range_x) x0(0) s(2) nose
    matrix RANGE_FRESH=e(diagnostics)
    // w_i=(12-i)/66. Absolute contribution coefficients are65 at row11,
    // and11,10,...,1,0 elsewhere; total130 and largest five sum103.
    test_assert "`family': fresh exceptional contributions have exact top-five share103/130" RANGE_FRESH[1,24]==1 & abs(RANGE_FRESH[1,25]-103/130)<1e-15
    tempfile range_model
    quietly estimates save `range_model', replace
    foreach state in fresh cold {
        if "`state'"=="cold" {
            preserve
            clear
            set obs 3
            gen byte unrelated=_n
            label data "Independent caller data"
            mata: mata clear
            quietly estimates use `range_model'
        }
        tempfile range_data
        quietly save "`range_data'", replace
        _pmpe_snapshot
        quietly rdistnn_results, type(influence)
        matrix RANGE_INFLUENCE=r(data)
        // Exact finite weighted residual: (1/66)*(8e307+191/33*1e307).
        scalar RANGE_EXPECTED=7e307*(65/2178)
        test_assert "`family' `state': numeric influence retains finite dominant contribution" RANGE_INFLUENCE[11,7]<. & abs(RANGE_INFLUENCE[11,7]/RANGE_EXPECTED-1)<1e-15
        quietly `family'_dnn_plot, type(influence) nlabel(0) name(range_influence) replace nodraw
        matrix RANGE_PLOT=r(plot_data)
        test_assert "`family' `state': influence graph retains finite dominant point" RANGE_PLOT[11,7]<. & abs(RANGE_PLOT[11,7]/RANGE_EXPECTED-1)<1e-15
        _pmpe_unchanged, data("`range_data'")
        local unchanged=r(ok)
        test_assert "`family' `state': influence views preserve complete e() and caller data" `unchanged'
        if "`state'"=="cold" restore
    }
    quietly rdistnn_predict, x0(0) post
    matrix RANGE_RESTORED=e(diagnostics)
    test_assert "`family': restored exceptional contributions have the same exact share" RANGE_RESTORED[1,24]==1 & abs(RANGE_RESTORED[1,25]-103/130)<1e-15
    quietly `family'_dnn range_signal, `extra' xcols(range_x) x0(0) s(12) nose
    matrix RANGE_ZERO=e(diagnostics)
    test_assert "`family': inactive extreme rows leave zero concentration available" RANGE_ZERO[1,24]==1 & missing(RANGE_ZERO[1,25])
    quietly rdistnn_predict, x0(0) post
    matrix RANGE_ZERO=e(diagnostics)
    test_assert "`family': restored zero concentration excludes inactive extremes" RANGE_ZERO[1,24]==1 & missing(RANGE_ZERO[1,25])
    clear
    set obs 6
    gen double range_x=_n
    gen double range_signal=cond(mod(_n,2),(-8e307)*2,8e307)
    gen double range_den=1
    quietly `family'_dnn range_signal, `extra' xcols(range_x) x0(0) s(1) nose
    matrix RANGE_TOTAL=e(diagnostics)
    // Each true absolute contribution is2e307: total1.2e308 exceeds Stata's
    // positive range, while concentration is exactly five sixths.
    test_assert "`family': overflowing contribution total retains exact5/6 share" RANGE_TOTAL[1,24]==1 & abs(RANGE_TOTAL[1,25]-5/6)<1e-15
    quietly rdistnn_predict, x0(0) post
    matrix RANGE_TOTAL=e(diagnostics)
    test_assert "`family': restored overflowing total retains exact5/6 share" RANGE_TOTAL[1,24]==1 & abs(RANGE_TOTAL[1,25]-5/6)<1e-15
}
test_section "tiny-weight recovered concentration retains decisive mass"
capture mata: RANGE_TINY=rdistnn_score_share((rdistnn_power2(-1074)\J(6,1,1)),((-8e307)*2\J(6,1,rdistnn_power2(-50))),(8e307\J(6,1,0))); st_numscalar("RANGE_TINY_AVAILABLE",RANGE_TINY[1]); st_numscalar("RANGE_TINY_SHARE",RANGE_TINY[2])
local range_rc=_rc
test_assert "tiny-weight recovered mass matches independent Fraction ratio" `range_rc'==0 & RANGE_TINY_AVAILABLE==1 & abs(RANGE_TINY_SHARE-0.7273363439834166)<1e-15

mata {
    RANGE_INACTIVE=rdistnn_score_share((0\1),( .\1),( .\0))
    RANGE_PRODUCT=rdistnn_score_share((1\0),(0\0),(8e307*8e307\0))
    st_numscalar("RANGE_INACTIVE_OK",RANGE_INACTIVE[1]==1 & RANGE_INACTIVE[2]==1)
    st_numscalar("RANGE_PRODUCT_OK",RANGE_PRODUCT[1]==0 & RANGE_PRODUCT[2]>=.)
}
mata {
    RANGE_SIGNED=J(0,1,.)
    RANGE_REUSE=rdistnn_score_share((rdistnn_power2(-1074)\J(6,1,1)),((-8e307)*2\J(6,1,rdistnn_power2(-50))),(8e307\J(6,1,0)),RANGE_SIGNED)
    st_numscalar("RANGE_SIGNED_FIRST",RANGE_SIGNED[1])
}
test_assert "tiny recovered dimensional contribution matches independent Fraction value" abs(RANGE_SIGNED_FIRST/-1.1857575500189916e-15-1)<1e-15

test_assert "inactive nonfinite operands do not enter concentration" RANGE_INACTIVE_OK
test_assert "genuinely nonfinite rounded paired products remain unavailable" RANGE_PRODUCT_OK

_rdistnn_test_e_snapshot_drop, name(PMPE_ESTATE)
macro drop PMPE_*
