// Observation/subset order, equal-distance weights and serialized predictions.
// Input index orders ranks within ties; it does not determine localization mass.
version 17
do "tests/test_framework.do"
test_section "exact-distance averaging with observation-order ranks"
clear
input double tx1 double tx2 double ty
1 0 40
-1 0 10
0 1 30
0 -1 20
end
gen long tid=_n
gen double tone=1
gen double tzero=0
gen double td=2*ty
matrix TT=(0,0\0,0)
matrix TR=(1\2\3\4)
foreach engine in mata plugin {
    global rdistnn_distance_engine `engine'
    foreach metric in l1 l2 l1_scaled l2_scaled mahalanobis {
        local extra ""
        if inlist("`metric'","l1_scaled","l2_scaled") local extra "scalerule(custom) scalevec(1 1)"
        quietly weights_dnn tx1 tx2, x0matrix(TT) s(2) distance(`metric') `extra'
        matrix TW=r(weights)
        matrix TRANK=r(ranks)
        matrix TID=(1\2\3\4)
        matrix TE=r(neff_targets)
        mata: st_numscalar("TIE_OK",all(st_matrix("TRANK"):==(st_matrix("TID"),st_matrix("TID"))))
        test_assert "`engine' `metric': exact ties rank by original observation IDs" TIE_OK
        test_assert "`engine' `metric': averaging expands support and changes actual ESS" abs(TW[1,1]-.25)<1e-15 & abs(TW[2,1]-.25)<1e-15 & abs(TW[3,1]-.25)<1e-15 & abs(TW[4,1]-.25)<1e-15 & abs(TE[1,1]-4)<1e-14
        test_assert "`engine' `metric': requested distance engine ran" "$rdistnn_last_distance_engine"=="`engine'"
    }
}
macro drop rdistnn_distance_engine
mata {
    TX=st_data(.,"tx1 tx2")
    // Reflection changes projected coordinate order, but every distance ties.
    TG=(1,2\0,0\-1,0\0,-1)
    st_numscalar("TIE_RAW",all(rdistnn_localization_ranks(rdistnn_distance_matrix(TX,(0,0),"mahalanobis",TG)):==st_matrix("TR")))
    st_numscalar("TIE_SUB",all(rdistnn_localization_ranks(rdistnn_distance_matrix(TX[(1\2\3),],(0,0))):==(1::3)))
    st_numscalar("TIE_DUP",all(rdistnn_localization_ranks(rdistnn_distance_matrix(TX\TX[1,],(0,0))):==(1::5)))
    st_numscalar("TIE_EXPLICIT",all(rdistnn_weights_from_ranks((1::4),1,4,"exact"):==(1\0\0\0)))
}
test_assert "Mahalanobis ties use original rows after whitening reflections" TIE_RAW
test_assert "subsets use original row order within exact-distance groups" TIE_SUB
test_assert "distinct and identical tied profiles use the same row order" TIE_DUP
test_assert "explicit supplied ranks remain authoritative" TIE_EXPLICIT

foreach family in np cate clate plqr tdnn {
    local cmd "`family'_dnn"
    local input "ty"
    local extra "s(3)"
    local targets "x0matrix(TT)"
    local expected=25
    if "`family'"=="cate" local extra "`extra' pseudooutcome"
    if "`family'"=="clate" local extra "`extra' pseudooutcome denominator(tone)"
    if "`family'"=="plqr" {
        local input "tone"
        local extra "s(3) d(td) md(tzero) g(tzero) thetaeval(0)"
        local targets "x0(0 0)"
    }
    if "`family'"=="tdnn" {
        local cmd "np_tdnn"
        local extra "s1(1) s2(3)"
        local expected=25
    }
    local xopt "xcols(tx1 tx2)"
    if "`family'"=="plqr" local xopt "x(tx1 tx2)"
    quietly `cmd' `input', `xopt' `targets' `extra'
    matrix TB=e(estimates)
    matrix TV=e(estimate_vcov)
    test_assert "`family': all-equal distinct distances give sample-mean jackknife variance" abs(TV[1,1]-125/3)<1e-11
    test_assert "`family': fresh point estimate follows canonical ties" abs(TB[1,1]-`expected')<1e-12
    test_assert "`family': fresh fit records canonical tie provenance" "`e(tie_rule)'"=="exact_distance_block_averaging_v1"
    tempfile cold
    quietly estimates save "`cold'", replace
    preserve
    gsort -tid
    quietly `cmd' `input', `xopt' `targets' `extra'
    matrix TP=e(estimates)
    test_assert "`family': reordering distinct covariates preserves the tied fit" abs(TP[1,1]-TB[1,1])<1e-12
    restore
    preserve
    clear
    mata: mata clear
    quietly estimates use "`cold'"
    quietly rdistnn_predict, x0matrix(TT) sefit post
    matrix TP=e(estimates)
    test_assert "`family': cold saved prediction keeps canonical ties and target identities" abs(TP[1,1]-`expected')<1e-12 & TP[1,1]==TP[1,2] & "`e(tie_rule)'"=="exact_distance_block_averaging_v1"
    restore
}

foreach family in np cate clate {
    local extra ""
    if "`family'"=="cate" local extra "pseudooutcome"
    if "`family'"=="clate" local extra "pseudooutcome denominator(tone)"
    quietly `family'_dnn ty, xcols(tx1 tx2) x0(0 0) `extra' tuningmode(grid) sgrid(2 3 2) sdefault(3)
    quietly rdistnn_candidate, candidate(1) post
    matrix TB=e(estimates)
    test_assert "`family': candidate extraction reuses canonical ranks without changing ESS" abs(TB[1,1]-25)<1e-12 & abs(e(neff)-4)<1e-14 & "`e(tie_rule)'"=="exact_distance_block_averaging_v1"
}

quietly np_dnn ty, xcols(tx1 tx2) x0(0 0) s(3) nose
quietly rdistnn_predict, x0(0 .5) post
matrix TN=e(estimates)
test_assert "new targets recompute metric distance groups from retained covariates" abs(TN[1,1]-28.75)<1e-12
