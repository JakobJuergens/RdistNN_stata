version 17
do "tests/helpers/fitted-data.do"
do "tests/helpers/e-state.do"
test_section "shared complete e() comparator owns estimation state"
quietly np_dnn fr_y, xcols(fr_x) x0(0) s(3)
estimates store STATE_REFERENCE
tempfile serialized
quietly estimates save "`serialized'", replace
quietly estimates use "`serialized'"
estimates store STATE_CURRENT
_rdistnn_test_e_compare STATE_REFERENCE, numeric
local equal=r(equal)
test_assert "numeric e() comparison ignores field enumeration order after serialization" `equal'
_rdistnn_test_e_compare STATE_CURRENT
local equal=r(equal)
test_assert "complete e() comparison preserves the active result" `equal'
capture program drop _state_test_mutate
program define _state_test_mutate, eclass
    syntax, KIND(string)
    if "`kind'"=="macro" ereturn local extra "descriptive metadata"
    if "`kind'"=="scalar" ereturn scalar extra=42
    if "`kind'"=="sample" {
        tempvar sample
        quietly gen byte `sample'=e(sample) & _n>1
        ereturn repost, esample(`sample')
    }
    if "`kind'"=="stripe" {
        tempname b
        matrix `b'=e(estimates)
        matrix colnames `b'=different_stripe
        ereturn matrix estimates=`b'
    }
end
foreach kind in macro scalar sample stripe {
    quietly estimates restore STATE_REFERENCE
    _state_test_mutate, kind(`kind')
    _rdistnn_test_e_compare STATE_REFERENCE
    local equal=r(equal)
    test_assert "complete e() comparison detects changed `kind'" !`equal'
    estimates store STATE_MUTATED
    _rdistnn_test_e_compare STATE_MUTATED
    local equal=r(equal)
    test_assert "comparison preserves changed `kind' state on return" `equal'
    if "`kind'"!="stripe" {
        _rdistnn_test_e_compare STATE_REFERENCE, numeric
        local equal=r(equal)
        test_assert "numeric mode compares matrices independently of descriptive `kind'" `equal'
    }
    estimates drop STATE_MUTATED
}
quietly estimates restore STATE_REFERENCE
gen str8 caller_text=cond(mod(_n,2),"odd","even")
gen byte caller_class=mod(_n,2)
label data "Complete e helper caller data"
label variable fr_x "caller covariate"
format fr_y %12.4f
label define STATE_CLASS 0 "control" 1 "treated"
label values caller_class STATE_CLASS
do "tests/helpers/caller-data.do"
_rdistnn_test_data_snapshot
tempfile caller
quietly save "`caller'", replace
_rdistnn_test_e_compare STATE_REFERENCE
capture cf _all using "`caller'"
local rc=_rc
test_assert "complete e() comparator leaves caller data and sample marker ownership intact" `rc'==0
_rdistnn_test_e_snapshot, name(STATE_INERT)
capture cf _all using "`caller'"
local rc=_rc
test_assert "inert complete e() snapshot creates no dataset marker variables" `rc'==0
_rdistnn_test_data_compare, data("`caller'")
local equal=r(equal)
test_assert "inert capture preserves caller values, variable order, storage, labels and formats" `equal'
_rdistnn_test_e_compare STATE_REFERENCE
local equal=r(equal)
test_assert "inert capture preserves every active e() field and e(sample)" `equal'
_rdistnn_test_e_compare, snapshot(STATE_INERT)
local equal=r(equal)
test_assert "inert snapshot compares every e() field and the sample exactly" `equal'
capture noisily _rdistnn_test_e_compare STATE_NO_SUCH_RESULT
local rc=_rc
test_assert "missing named reference returns a failure" `rc'!=0
_rdistnn_test_e_compare, snapshot(STATE_INERT)
local equal=r(equal)
_rdistnn_test_data_compare, data("`caller'")
local same_data=r(equal)
test_assert "failed named comparison preserves complete caller e() and data" `equal' & `same_data'
_rdistnn_test_e_snapshot_drop, name(STATE_INERT)
estimates drop STATE_CURRENT STATE_REFERENCE
program drop _state_test_mutate

// Caller Mata bindings and absence are state too. Keep the snapshot local to
// this test function while exercising real public commands, including errors.
capture mata: mata drop rdistnn_test_scratch_scope()
mata:
void rdistnn_test_scratch_scope(string scalar command, string rowvector names,
    real scalar present, string scalar output)
{
    transmorphic scalar saved
    real rowvector existed
    real scalar i,rc,same,no_new
    pointer() scalar p
    string colvector before,after
    saved=asarray_create("string",1); existed=J(1,cols(names),0)
    for(i=1;i<=cols(names);i++) {
        p=findexternal(names[i])
        if(p!=NULL) {
            existed[i]=1
            asarray(saved,names[i],(*p))
        }
        rmexternal(names[i])
        if(present) {
            p=crexternal(names[i])
            *p=(42,i)
        }
    }
    before=sort(direxternal("*"),1)
    rc=_stata(command,1)
    same=1
    for(i=1;i<=cols(names);i++) {
        p=findexternal(names[i])
        if(present) {
            if(p==NULL) same=0
            else if(rows((*p))!=1 | cols((*p))!=2) same=0
            else if(any((*p):!=(42,i))) same=0
        }
        else if(p!=NULL) same=0
    }
    after=sort(direxternal("*"),1)
    no_new=rows(before)==rows(after)
    if(no_new & rows(before)>0) no_new=all(before:==after)
    for(i=1;i<=cols(names);i++) {
        rmexternal(names[i])
        if(existed[i]) {
            p=crexternal(names[i])
            *p=asarray(saved,names[i])
        }
    }
    st_matrix(output,(rc,same,no_new))
}
end

test_section "public numerical and plot work leaves Mata scratch bindings local"
clear
set obs 12
gen double scratch_x=mod(_n,6)-3
gen double scratch_y=_n^2
gen double scratch_den=1+_n/20
gen double scratch_zero=0
matrix SCRATCH_TARGETS=(0\1)
local scratch_names "_cl_A _cl_CV _cl_DX _cl_V _cl_W _cl_X _cl_X0 _cl_b _cl_c _cl_e _cl_g _cl_geometry _cl_grid _cl_grid_W _cl_grid_diag _cl_gridd _cl_gridse _cl_gridstatus _cl_gridy _cl_groups _cl_h _cl_keep _cl_meta _cl_prior _cl_ranks _cl_signals _cl_status _ct_DX _ct_E _ct_T _ct_V _ct_Vg _ct_W _ct_Wg _ct_X _ct_X0 _ct_candidates _ct_g _ct_geometry _ct_grid_diag _ct_grid_est _ct_grid_se _ct_groups _ct_h _ct_neff _ct_prior _ct_ranks _ct_rho _ct_sg _ct_target_diag _ct_tau _np_DX _np_V _np_Vg _np_W _np_Wg _np_X _np_X0 _np_candidates _np_g _np_geometry _np_grid_diag _np_grid_est _np_grid_se _np_groups _np_h _np_mu _np_neff _np_prior _np_ranks _np_sg _np_target_diag _np_y _pl_D _pl_DX _pl_V _pl_X _pl_below _pl_g0 _pl_geometry _pl_groups _pl_md _pl_moment _pl_neff _pl_ranks _pl_resid _pl_score _pl_theta _pl_w _pl_x0 _pl_y _rdnn_av _rdnn_ids _rdnn_pp _rdnn_ss _td_DX _td_R _td_W _td_X _td_X0 _td_geom _td_groups _td_plot _td_table _td_y _wd_D _wd_G _wd_R _wd_W _wd_X _wd_X0 _wd_groups"
tempname scratch_result
foreach present in 0 1 {
    foreach family in np cate clate plqr tdnn weights tdweights {
        local command "np_dnn scratch_y, xcols(scratch_x) x0matrix(SCRATCH_TARGETS) s(4)"
        if "`family'"=="cate" local command "cate_dnn scratch_y, pseudooutcome xcols(scratch_x) x0matrix(SCRATCH_TARGETS) s(4)"
        if "`family'"=="clate" local command "clate_dnn scratch_y, pseudooutcome denominator(scratch_den) xcols(scratch_x) x0matrix(SCRATCH_TARGETS) tuningmode(grid) sgrid(2 4) sdefault(4)"
        if "`family'"=="plqr" local command "plqr_dnn scratch_y, x(scratch_x) x0(0) d(scratch_den) g(scratch_zero) md(scratch_zero) thetaeval(0) s(4)"
        if "`family'"=="tdnn" local command "np_tdnn scratch_y, xcols(scratch_x) x0matrix(SCRATCH_TARGETS) s1(2) s2(4)"
        if "`family'"=="weights" local command "weights_dnn scratch_x, x0matrix(SCRATCH_TARGETS) s(4)"
        if "`family'"=="tdweights" local command "weights_tdnn scratch_x, x0matrix(SCRATCH_TARGETS) s1(2) s2(4)"
        foreach failure in 0 1 {
            local call "quietly `command'"
            local expected 0
            if `failure' {
                local call "`call' distance(l2_scaled) scalerule(custom) scalevec(1 2)"
                local expected 198
            }
            mata: rdistnn_test_scratch_scope(st_local("call"),tokens(st_local("scratch_names")),`present',st_local("scratch_result"))
            test_assert "`family'/`present'/`failure': public command preserves its return code" `scratch_result'[1,1]==`expected'
            test_assert "`family'/`present'/`failure': existing or absent scratch bindings survive" `scratch_result'[1,2]==1
            test_assert "`family'/`present'/`failure': no new Mata arrays remain" `scratch_result'[1,3]==1
        }
    }
    foreach plot in tdnn np clate {
        if "`plot'"=="tdnn" quietly np_tdnn scratch_y, xcols(scratch_x) x0matrix(SCRATCH_TARGETS) s1(2) s2(4)
        if "`plot'"=="np" quietly np_dnn scratch_y, xcols(scratch_x) x0matrix(SCRATCH_TARGETS) s(4)
        if "`plot'"=="clate" quietly clate_dnn scratch_y, pseudooutcome denominator(scratch_den) xcols(scratch_x) x0matrix(SCRATCH_TARGETS) tuningmode(grid) sgrid(2 4) sdefault(4)
        local command "np_dnn_plot, name(scratch_scope) nodraw"
        if "`plot'"=="tdnn" local command "np_tdnn_plot, type(weights) name(scratch_scope) nodraw"
        if "`plot'"=="clate" local command "clate_dnn_plot, candidates(1 2) name(scratch_scope) nodraw"
        capture graph drop scratch_scope
        foreach failure in 0 1 {
            local call "quietly `command'"
            local expected=cond(`failure',110,0)
            mata: rdistnn_test_scratch_scope(st_local("call"),tokens(st_local("scratch_names")),`present',st_local("scratch_result"))
            test_assert "`plot'/`present'/`failure': plotting preserves its return code" `scratch_result'[1,1]==`expected'
            test_assert "`plot'/`present'/`failure': plotting preserves scratch bindings" `scratch_result'[1,2]==1
            test_assert "`plot'/`present'/`failure': plotting retains no Mata arrays" `scratch_result'[1,3]==1
        }
        graph drop scratch_scope
    }
}
mata: mata drop rdistnn_test_scratch_scope()
