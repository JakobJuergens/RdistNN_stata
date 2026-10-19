// Actual Rcpp/R references, native bridge validation and forced route evidence.
version 17
quietly do "src/_rdistnn_mata.ado"
do "tests/test_framework.do"
quietly do "tests/fixtures/r_native_distance_fixtures.do"
test_section "native distances: versioned protocol and actual R references"
matrix ND_API=(73)
capture plugin call rdistnn_plugin, distance_api ND_API
local api_rc=_rc
test_assert "distance API version1 is actually callable" `api_rc'==0 & ND_API[1,1]==1
// Old plugins must visibly fail this test rather than skip the new contract.
if `api_rc'!=0 exit 0
capture program drop _nd_case
program define _nd_case
    args id metric rule xx target
    local inv "none"
    local distance "`metric'"
    if "`rule'"!="none" {
        local inv "ND_`id'_inv"
        local distance "`metric'_scaled"
    }
    matrix ND_OUT=J(rowsof(`xx'),rowsof(`target'),73)
    capture plugin call rdistnn_plugin, distance_batch `metric' `xx' `target' `inv' ND_OUT
    local native_rc=_rc
    test_assert "`id': native distance_batch succeeds" `native_rc'==0
    if `native_rc'!=0 exit
    mata {
        nd_c=st_matrix("ND_`id'_C"); nd_r=st_matrix("ND_`id'_R")
        st_numscalar("nd_err",max(abs(st_matrix("ND_OUT")-nd_c):/(1:+abs(nd_c))))
        st_numscalar("nd_rank",max(abs(rdistnn_ranks_from_dmat(st_matrix("ND_OUT"))-st_matrix("ND_`id'_rankC"))))
        st_numscalar("nd_rdiff",sum(st_matrix("ND_`id'_rankC"):!=st_matrix("ND_`id'_rankR")))
    }
    test_assert "`id': native distances match actual Rcpp" scalar(nd_err)<1e-14
    test_assert "`id': native stable ranks match actual Rcpp" scalar(nd_rank)==0
    test_assert "`id': R native/fallback rank discrepancies are retained" scalar(nd_rdiff)==ND_`id'_rankdiff[1,1]
    test_assert "`id': repaired R fallback preserves actual native ranks" scalar(nd_rdiff)==0
    // Ordinary fallback distances agree numerically, but never substitute its
    // potentially different ranks for the Rcpp arithmetic contract.
    if strpos("`id'","offset")==0 {
        mata: st_numscalar("nd_err",max(abs(nd_c-nd_r):/(1:+abs(nd_c))))
        test_assert "`id': maintained R fallback numeric tolerance" scalar(nd_err)<1e-14
    }
    foreach engine in mata plugin {
        global rdistnn_distance_engine `engine'
        mata: nd_G=("`rule'"=="none" ? J(0,0,.) : st_matrix("ND_`id'_scale")); nd_D=rdistnn_distance_matrix(st_matrix("`xx'"),st_matrix("`target'"),"`distance'",nd_G)
        local observed "$rdistnn_last_distance_engine"
        mata {
            st_numscalar("nd_err",max(abs(nd_D-nd_c):/(1:+abs(nd_c))))
            st_numscalar("nd_rank",max(abs(rdistnn_ranks_from_dmat(nd_D)-st_matrix("ND_`id'_rankC"))))
            nd_perm=(rows(st_matrix("`target'"))..1)
            nd_P=rdistnn_distance_matrix(st_matrix("`xx'"),st_matrix("`target'")[nd_perm,],"`distance'",nd_G)
            st_numscalar("nd_perm",max(abs(nd_P-nd_D[,nd_perm])))
            nd_S=rdistnn_distance_matrix(st_matrix("`xx'"),st_matrix("`target'")[1,],"`distance'",nd_G)
            st_numscalar("nd_scalar",max(abs(nd_S-nd_D[,1])))
        }
        test_assert "`id'/`engine': dispatcher observed intended route" "`observed'"=="`engine'"
        test_assert "`id'/`engine': actual Rcpp distances" scalar(nd_err)<1e-14
        test_assert "`id'/`engine': actual Rcpp stable ranks" scalar(nd_rank)==0
        test_assert "`id'/`engine': reordered and duplicated targets retain columns" scalar(nd_perm)==0
        test_assert "`id'/`engine': scalar equals batch" scalar(nd_scalar)==0
    }
end
foreach metric in l1 l2 {
    _nd_case `metric' `metric' none ND_X ND_X0
    foreach rule in sd mad iqr range custom {
        mata: st_numscalar("nd_scale",mreldif(rdistnn_fit_geometry(st_matrix("ND_X"),"`metric'_scaled","`rule'",("`rule'"=="custom" ? (-2,0,4) : J(1,0,.)),1e-8),st_matrix("ND_`metric'_`rule'_scale")))
        test_assert "`metric'/`rule': fitted scales match sibling R" scalar(nd_scale)<1e-14
        _nd_case `metric'_`rule' `metric' `rule' ND_X ND_X0
    }
    _nd_case `metric'_tie `metric' none ND_TX ND_T0
    _nd_case `metric'_offset `metric' custom ND_OX ND_O0
    _nd_case `metric'_extreme `metric' custom ND_EX ND_E0
}
_nd_case l1_accum l1 none ND_AX1 ND_A0
_nd_case l2_accum l2 none ND_AX2 ND_A0
_nd_case l2_rounding l2 none ND_FX ND_F0
test_assert "ordinary scaled R rank discrepancy is repaired" ND_l2_custom_rankdiff[1,1]==0
test_assert "large-offset scaled R rank discrepancy is repaired" ND_l2_offset_rankdiff[1,1]==0

test_section "native distances: malformed bridge inputs preserve all output"
matrix ND_BAD0=(0,0)
matrix ND_BADX=ND_X
matrix ND_BADX[2,1]=.
matrix ND_BADT=ND_X0
matrix ND_BADT[2,2]=.a
matrix ND_BADINV=(1,.,1)
matrix ND_SHORT=(1,1)
matrix ND_COL=(1\1\1)
matrix ND_ZERO=(1,0,1)
matrix ND_NEGATIVE=(1,-1,1)
foreach bad in metric width xmissing tmissing invmissing invshort invcolumn invzero invnegative output missing extra nonexistent {
    matrix ND_OUT=J(6,4,73)
    local op "distance_batch l2 ND_X ND_X0 none ND_OUT"
    if "`bad'"=="metric" local op "distance_batch l3 ND_X ND_X0 none ND_OUT"
    if "`bad'"=="width" local op "distance_batch l2 ND_X ND_BAD0 none ND_OUT"
    if "`bad'"=="xmissing" local op "distance_batch l2 ND_BADX ND_X0 none ND_OUT"
    if "`bad'"=="tmissing" local op "distance_batch l2 ND_X ND_BADT none ND_OUT"
    if "`bad'"=="invmissing" local op "distance_batch l2 ND_X ND_X0 ND_BADINV ND_OUT"
    if "`bad'"=="invshort" local op "distance_batch l2 ND_X ND_X0 ND_SHORT ND_OUT"
    if "`bad'"=="invcolumn" local op "distance_batch l2 ND_X ND_X0 ND_COL ND_OUT"
    if "`bad'"=="invzero" local op "distance_batch l2 ND_X ND_X0 ND_ZERO ND_OUT"
    if "`bad'"=="invnegative" local op "distance_batch l2 ND_X ND_X0 ND_NEGATIVE ND_OUT"
    if "`bad'"=="output" {
        matrix ND_OUT=J(6,3,73)
    }
    if "`bad'"=="missing" local op "distance_batch l2 ND_X ND_X0 ND_OUT"
    if "`bad'"=="extra" local op "distance_batch l2 ND_X ND_X0 none ND_OUT surplus"
    if "`bad'"=="nonexistent" local op "distance_batch l2 ND_absent ND_X0 none ND_OUT"
    capture plugin call rdistnn_plugin, `op'
    local bad_rc=_rc
    mata: st_numscalar("nd_untouched",all(st_matrix("ND_OUT"):==73))
    test_assert "`bad': malformed native input fails before output writes" `bad_rc'!=0 & scalar(nd_untouched)
}
matrix ND_HUGE=(1,1 \ 8e307,8e307)
matrix ND_ORIGIN=(0,0)
matrix ND_OUT=J(2,1,73)
capture plugin call rdistnn_plugin, distance_batch l2 ND_HUGE ND_ORIGIN none ND_OUT
local overflow_rc=_rc
mata: st_numscalar("nd_untouched",all(st_matrix("ND_OUT"):==73))
test_assert "overflow native distance fails with untouched output" `overflow_rc'!=0 & scalar(nd_untouched)
matrix ND_API=J(1,2,73)
capture plugin call rdistnn_plugin, distance_api ND_API
local bad_rc=_rc
mata: st_numscalar("nd_untouched",all(st_matrix("ND_API"):==73))
test_assert "malformed capability output preserved" `bad_rc'!=0 & scalar(nd_untouched)

test_section "native distances: matrices above matsize and geometry boundaries"
local nd_matsize=c(matsize)
mata {
    nd_n=st_numscalar("c(matsize)")+1
    nd_tall=((1::nd_n),sin(1::nd_n))
    nd_wide=((1::nd_n),cos(1::nd_n))
    nd_small=(0,0\1,1\2,2)
}
foreach metric in l1 l2 {
    global rdistnn_distance_engine mata
    mata: nd_tref=rdistnn_distance_matrix(nd_tall,nd_small,"`metric'");nd_wref=rdistnn_distance_matrix(nd_small,nd_wide,"`metric'")
    global rdistnn_distance_engine plugin
    mata {
        nd_large=rdistnn_distance_matrix(nd_tall,nd_small,"`metric'")
        st_numscalar("nd_tallerr",max(abs(nd_large-nd_tref)))
        nd_large=rdistnn_distance_matrix(nd_small,nd_wide,"`metric'")
        st_numscalar("nd_wideerr",max(abs(nd_large-nd_wref)))
    }
    test_assert "`metric': native tall and wide bridges exceed matsize" scalar(nd_tallerr)<1e-12 & scalar(nd_wideerr)<1e-12
}
test_assert "native distance bridge preserves user matsize" c(matsize)==`nd_matsize'
foreach engine in mata plugin {
    global rdistnn_distance_engine `engine'
    foreach bad in short zero missing {
        local geo "(1,2)"
        if "`bad'"=="zero" local geo "(1,0,1)"
        if "`bad'"=="missing" local geo "(1,.,1)"
        capture mata: nd_bad=rdistnn_distance_matrix(st_matrix("ND_X"),st_matrix("ND_X0"),"l2_scaled",`geo')
        test_assert "`engine'/`bad': invalid frozen geometry fails" _rc!=0
    }
    capture mata: nd_bad=rdistnn_distance_matrix(J(0,3,.),st_matrix("ND_X0"),"l2")
    test_assert "`engine': empty training matrix rejected" _rc!=0
    capture mata: nd_bad=rdistnn_distance_matrix(st_matrix("ND_X"),J(0,3,.),"l2")
    test_assert "`engine': empty target matrix rejected" _rc!=0
}

test_section "native distances: dispatch controls and sample mapping"
global rdistnn_distance_engine auto
capture macro drop rdistnn_use_plugin
mata: nd_D=rdistnn_distance_matrix(st_matrix("ND_X"),st_matrix("ND_X0"),"l2")
test_assert "auto avoids native bridge overhead for small calls" "$rdistnn_last_distance_engine"=="mata"
global rdistnn_use_plugin 0
mata: nd_D=rdistnn_distance_matrix(st_matrix("ND_X"),st_matrix("ND_X0"),"l2")
test_assert "auto respects use_plugin0" "$rdistnn_last_distance_engine"=="mata"
global rdistnn_use_plugin 1
mata: nd_D=rdistnn_distance_matrix(st_matrix("ND_X"),st_matrix("ND_X0"),"l2")
test_assert "auto uses supported native capability" "$rdistnn_last_distance_engine"=="plugin"
capture macro drop rdistnn_use_plugin
mata {
    nd_batch_X=J(5000,16,0)
    nd_batch_T=J(256,16,1)
    nd_D=rdistnn_distance_matrix(nd_batch_X,nd_batch_T,"l2")
    st_numscalar("nd_batch_ok",all(nd_D:==4))
}
test_assert "auto uses native for a substantial batch" "$rdistnn_last_distance_engine"=="plugin" & scalar(nd_batch_ok)
global rdistnn_use_plugin 1
global rdistnn_distance_engine mata
mata: nd_D=rdistnn_distance_matrix(st_matrix("ND_X"),st_matrix("ND_X0"),"l2")
test_assert "explicit Mata overrides use_plugin1" "$rdistnn_last_distance_engine"=="mata"
global rdistnn_distance_engine invalid
capture mata: nd_D=rdistnn_distance_matrix(st_matrix("ND_X"),st_matrix("ND_X0"),"l2")
test_assert "invalid engine fails early" _rc==198
global rdistnn_distance_engine auto
capture macro drop rdistnn_use_plugin
local original_ok "$rdistnn_plugin_ok"
global rdistnn_plugin_ok 0
mata: nd_D=rdistnn_distance_matrix(st_matrix("ND_X"),st_matrix("ND_X0"),"l2")
test_assert "missing plugin auto falls back" "$rdistnn_last_distance_engine"=="mata"
global rdistnn_use_plugin 1
capture mata: nd_D=rdistnn_distance_matrix(st_matrix("ND_X"),st_matrix("ND_X0"),"l2")
test_assert "missing plugin auto with shared forcing fails" _rc!=0
capture macro drop rdistnn_use_plugin
global rdistnn_distance_engine plugin
capture mata: nd_D=rdistnn_distance_matrix(st_matrix("ND_X"),st_matrix("ND_X0"),"l2")
test_assert "missing plugin forcing fails explicitly" _rc!=0
global rdistnn_plugin_ok `original_ok'
global rdistnn_distance_engine plugin
mata: nd_D=rdistnn_distance_matrix(st_matrix("ND_X"),st_matrix("ND_X0"),"mahalanobis")
test_assert "Mahalanobis remains Mata under forced native selector" "$rdistnn_last_distance_engine"=="mata"
clear
svmat double ND_X, names(ndx)
gen double nd_y=ndx1+2*ndx2
replace ndx1=. in 2
foreach metric in l1 l2 l1_scaled l2_scaled {
    global rdistnn_distance_engine mata
    quietly weights_dnn ndx1 ndx2 ndx3 if _n!=5, x0matrix(ND_X0) s(3) distance(`metric') scalerule(custom) scalevec(-2 0 4) scaleeps(1e-8)
    matrix ND_W=r(weights)
    matrix ND_RANK=r(ranks)
    matrix ND_IDS=r(row_ids)
    global rdistnn_distance_engine plugin
    quietly weights_dnn ndx1 ndx2 ndx3 if _n!=5, x0matrix(ND_X0) s(3) distance(`metric') scalerule(custom) scalevec(-2 0 4) scaleeps(1e-8)
    local observed "$rdistnn_last_distance_engine"
    mata: st_numscalar("nd_equal",max(abs(st_matrix("ND_W")-st_matrix("r(weights)")))+max(abs(st_matrix("ND_RANK")-st_matrix("r(ranks)")))+max(abs(st_matrix("ND_IDS")-st_matrix("r(row_ids)"))))
    test_assert "`metric': weights really use native route" "`observed'"=="plugin"
    test_assert "`metric': complete sample ranks/weights/row mappings preserved" scalar(nd_equal)==0 & ND_IDS[2,1]==3 & ND_IDS[4,1]==6
    capture weights_dnn ndx1 ndx2 ndx3, x0(0 2 5) s(3) distance(`metric') scalerule(custom) scalevec(1 2)
    test_assert "`metric': wrong custom scale dimension rejected" _rc==198
    capture weights_dnn ndx1 ndx2 ndx3, x0(0 2 5) s(3) distance(`metric') scaleeps(0)
    test_assert "`metric': nonpositive scale epsilon rejected" _rc==198
}
capture macro drop rdistnn_distance_engine
capture macro drop rdistnn_last_distance_engine
capture macro drop rdistnn_use_plugin
capture program drop _nd_case
