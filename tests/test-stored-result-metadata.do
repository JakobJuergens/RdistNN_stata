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

do "tests/helpers/fitted-data.do"
test_section "CLATE prediction covariance metadata follows requested policy"
quietly clate_dnn fr_y, denominator(fr_den) pseudooutcome xcols(fr_x) x0(0) s(4)
quietly rdistnn_predict, x0(1) post
test_assert "CLATE no-SE prediction clears all covariance metadata" "`e(se_method)'"=="none" & "`e(covariance)'"=="none" & "`e(jackknife)'"=="none"
quietly rdistnn_predict, x0(1) sefit post
test_assert "CLATE sefit enables covariance metadata on a no-SE parent" "`e(se_method)'"=="paired_rank_shift" & "`e(covariance)'"=="paired_frozen_score" & "`e(jackknife)'"=="fast"


test_section "candidate metadata uses current result ownership"
clear
set obs 12
gen double fx=(_n-6)/5
gen double fy=1+fx^2+.2*sin(_n)
gen double fd=.8+.05*fx
matrix FT=(-.4\ .4)
foreach family in np cate clate {
    local extra ""
    if "`family'"=="cate" local extra "pseudooutcome"
    if "`family'"=="clate" local extra "pseudooutcome denominator(fd)"
    local vce=cond("`family'"=="clate","paired_rank_shift","rank_shift")
    quietly `family'_dnn fy, xcols(fx) x0matrix(FT) tuningmode(grid) sgrid(2 3) sdefault(3) `extra'
    quietly rdistnn_candidate, candidate(1) post
    test_assert "`family': candidate retains the input variable name" "`e(depvar)'"=="fy"
    test_assert "`family': candidate retains the VCE method" "`e(vce)'"=="`vce'"
    quietly rdistnn_predict, x0(0) sefit post
    test_assert "`family': prediction retains the VCE method" "`e(vce)'"=="`vce'"
    quietly rdistnn_predict, x0(0) post
    test_assert "`family': point-only prediction records no VCE" "`e(vce)'"=="none"
    test_assert "`family': point-only prediction clears the VCE label" "`e(vcetype)'"==""
}

test_section "saved TDNN keeps signed concentration fields without probability ESS aliases"
clear
set obs 12
gen double sc_x=_n
gen double sc_y=sin(_n)
matrix SC_TARGETS=(0\13)
foreach targetset in scalar batch {
    local targets "x0(0)"
    if "`targetset'"=="batch" local targets "x0matrix(SC_TARGETS)"
    quietly np_tdnn sc_y, xcols(sc_x) `targets' s1(2) s2(5)
    matrix SC_DIAGNOSTICS=e(diagnostics)
    matrix SC_CONCENTRATION=e(neff_targets)
    foreach phase in fresh derived {
        if "`phase'"=="derived" quietly rdistnn_predict, `targets' sefit post
        local scalars : e(scalars)
        local matrices : e(matrices)
        foreach field in ess ess_min ess_max {
            local present : list field in scalars
            test_assert "TDNN `targetset'/`phase': omits probability-named `field'" !`present'
        }
        local field "ess_targets"
        local present : list field in matrices
        test_assert "TDNN `targetset'/`phase': omits probability-named weight ESS matrix" !`present'
        mata: st_numscalar("sc_unchanged",all(st_matrix("e(diagnostics)"):==st_matrix("SC_DIAGNOSTICS")) & max(abs(st_matrix("e(neff_targets)")-st_matrix("SC_CONCENTRATION")))<=8*epsilon(1)*max((1,max(abs(st_matrix("SC_CONCENTRATION"))))))
        test_assert "TDNN `targetset'/`phase': keeps signed diagnostics and concentration" scalar(sc_unchanged) & el(e(diagnostics),1,3)>0
    }
}
