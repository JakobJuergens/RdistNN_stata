// Corrupt schema-2 snapshots must fail before reuse or observation reconstruction.
version 17
do "tests/test_framework.do"
do "tests/helpers/e-state.do"
local smv_distance "$rdistnn_distance_engine"
local smv_jackknife "$rdistnn_jackknife_engine"
global rdistnn_distance_engine mata
global rdistnn_jackknife_engine mata
clear
set obs 8
gen double smv_x=_n-4.5
gen double smv_z=0
gen double smv_y=_n^2
gen double smv_den=1+_n/100
gen double smv_d=2+sin(_n)
gen double smv_d2=-3+cos(_n)
gen double smv_zero=0
matrix SMV_TARGETS=(0,0\2.5,0)
quietly np_dnn smv_y, xcols(smv_x smv_z) x0matrix(SMV_TARGETS) s(3)
estimates store SMV_PARENT
capture program drop _smv_corrupt
program define _smv_corrupt, eclass
    syntax, CASE(string)
    tempname value
    if inlist("`case'","negative_distance","missing_distance") {
        matrix `value'=e(model_distances)
        matrix `value'[1,1]=cond("`case'"=="negative_distance",-1,.)
        ereturn matrix model_distances=`value', copy
    }
    if inlist("`case'","rank_order","rank_tie","duplicate_rank") {
        matrix `value'=e(model_ranks)
        local a=1
        local b=2
        if "`case'"=="rank_tie" {
            local a=4
            local b=5
        }
        local rank=el(`value',`a',1)
        matrix `value'[`a',1]=el(`value',`b',1)
        if "`case'"!="duplicate_rank" matrix `value'[`b',1]=`rank'
        ereturn matrix model_ranks=`value', copy
    }
    if inlist("`case'","fractional_id","duplicate_id","missing_id","short_ids") {
        matrix `value'=e(model_row_ids)
        if "`case'"=="fractional_id" matrix `value'[1,1]=-1.5
        if "`case'"=="duplicate_id" matrix `value'[2,1]=el(`value',1,1)
        if "`case'"=="missing_id" matrix `value'[1,1]=.
        if "`case'"=="short_ids" matrix `value'=`value'[1..7,1]
        ereturn matrix model_row_ids=`value', copy
    }
    if "`case'"=="input_rows" {
        matrix `value'=e(model_inputs)
        matrix `value'=`value'[1..7,1...]
        ereturn matrix model_inputs=`value', copy
    }
    if "`case'"=="input_names" ereturn local model_vars "smv_x smv_z"
    if "`case'"=="signal_columns" {
        matrix `value'=e(model_signals), e(model_signals)
        ereturn matrix model_signals=`value', copy
    }
    if "`case'"=="signal_missing" {
        matrix `value'=e(model_signals)
        matrix `value'[1,1]=.
        ereturn matrix model_signals=`value', copy
    }
    if "`case'"=="method" ereturn local method "invalid"
    if "`case'"=="distance" ereturn local distance "invalid"
    if "`case'"=="count" ereturn scalar N=9
    if "`case'"=="dimension" ereturn scalar p=3
    if "`case'"=="target_count" ereturn scalar m=3
    if "`case'"=="legacy" ereturn scalar model_version=1
    if "`case'"=="geometry" {
        matrix `value'=(1,0)
        ereturn matrix model_geometry=`value', copy
    }
end

test_section "saved snapshot corruption: reject consistently and preserve caller state"
foreach case in negative_distance missing_distance rank_order rank_tie duplicate_rank ///
    fractional_id duplicate_id missing_id short_ids input_rows input_names ///
    signal_columns signal_missing method distance count dimension target_count legacy {
    quietly estimates restore SMV_PARENT
    _smv_corrupt, case(`case')
    estimates store SMV_CORRUPT_ACTIVE
    quietly estimates restore SMV_CORRUPT_ACTIVE
    _rdistnn_test_e_snapshot, name(SMV_CORRUPT)
    foreach route in cached new weights kernel profile {
        quietly estimates restore SMV_CORRUPT_ACTIVE
        if "`route'"=="cached" capture quietly rdistnn_predict, x0matrix(SMV_TARGETS) post
        if "`route'"=="new" capture quietly rdistnn_predict, x0(.25 0) post
        if "`route'"=="weights" capture quietly rdistnn_results, type(weights)
        if "`route'"=="kernel" capture quietly np_dnn_plot, type(weights) view(kernel) name(smv_invalid) replace nodraw
        if "`route'"=="profile" capture quietly np_dnn_plot, type(weights) name(smv_invalid) replace nodraw
        local rc=_rc
        test_assert "`case'/`route': corrupt snapshot is rejected" `rc'==498
        _rdistnn_test_e_compare, snapshot(SMV_CORRUPT)
        local same=r(equal)
        test_assert "`case'/`route': failure preserves complete e() and sample" `same' & _N==8 & smv_y[1]==1
    }
    estimates drop SMV_CORRUPT_ACTIVE
    _rdistnn_test_e_snapshot_drop, name(SMV_CORRUPT)
}

// Malformed frozen geometry is invalid even when the requested targets are cached.
quietly np_dnn smv_y, xcols(smv_x smv_z) x0matrix(SMV_TARGETS) s(3) distance(l2_scaled)
_smv_corrupt, case(geometry)
capture quietly rdistnn_predict, x0matrix(SMV_TARGETS)
local rc=_rc
test_assert "cached targets reject a nonpositive frozen scale" `rc'==498
capture quietly rdistnn_results, type(weights)
local rc=_rc
test_assert "weight extraction rejects a nonpositive frozen scale" `rc'==498

test_section "valid cold snapshots: all estimator families retain fitted geometry and weights"
replace smv_x=smv_x[1] in 2
foreach family in np tdnn cate clate plqr {
    local covariates "xcols(smv_x smv_z) x0matrix(SMV_TARGETS)"
    if "`family'"=="np" quietly np_dnn smv_y, `covariates' s(3) distance(l2_scaled)
    if "`family'"=="tdnn" quietly np_tdnn smv_y, `covariates' s1(2) s2(3) distance(l2_scaled)
    if "`family'"=="cate" quietly cate_dnn smv_y, pseudooutcome `covariates' s(3) distance(l2_scaled)
    if "`family'"=="clate" quietly clate_dnn smv_y, pseudooutcome denominator(smv_den) `covariates' s(3) distance(l2_scaled)
    if "`family'"=="plqr" quietly plqr_dnn smv_y, x(smv_x smv_z) x0(0 0) d(smv_d smv_d2) md(smv_zero smv_zero) g(smv_zero) thetaeval(0 0) s(3) distance(l2_scaled)
    matrix SMV_FIT=e(estimates)
    matrix SMV_GEOMETRY=e(model_geometry)
    matrix SMV_CURRENT=e(model_targets)
    tempfile cold
    quietly estimates save "`cold'", replace
    preserve
    clear
    mata: mata clear
    quietly estimates use "`cold'"
    quietly rdistnn_predict, x0matrix(SMV_CURRENT) sefit
    matrix SMV_GOT=r(fit)
    mata: st_numscalar("smv_same",all(vec(st_matrix("SMV_GOT"))':==st_matrix("SMV_FIT")) & all(st_matrix("SMV_GEOMETRY"):==st_matrix("e(model_geometry)")))
    test_assert "`family': cold cached prediction keeps estimates and frozen geometry" scalar(smv_same) & _N==0
    quietly rdistnn_results, type(weights)
    local rows=rowsof(r(data))
    test_assert "`family': valid cold observation extraction keeps every ID" `rows'==e(N)*e(m)
    quietly rdistnn_predict, x0(.25 0) sefit
    local columns=colsof(r(fit))
    test_assert "`family': valid cold prediction supports new targets" `columns'==1 & _N==0
    if inlist("`family'","cate","clate") {
        quietly rdistnn_results, type(influence)
        local rows=rowsof(r(data))
        test_assert "`family': valid cold influence retains every row" `rows'==e(N)*e(m)
        _smv_corrupt, case(negative_distance)
        capture quietly rdistnn_results, type(influence)
        local rc=_rc
        test_assert "`family': influence rejects corrupt saved distances" `rc'==498
    }
    if inlist("`family'","tdnn","plqr","clate") {
        quietly estimates use "`cold'"
        _smv_corrupt, case(signal_columns)
        capture quietly rdistnn_predict, x0matrix(SMV_CURRENT)
        local rc=_rc
        test_assert "`family': wrong signal component count is rejected" `rc'==498
    }
    restore
}
quietly estimates restore SMV_PARENT
estimates drop SMV_PARENT
program drop _smv_corrupt
global rdistnn_distance_engine "`smv_distance'"
global rdistnn_jackknife_engine "`smv_jackknife'"
