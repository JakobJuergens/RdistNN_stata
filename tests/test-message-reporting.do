version 17
do "tests/test_framework.do"
test_section "shared reporting channels and state"
clear
set obs 12
gen double x=_n
gen double y=2*x
quietly np_dnn y, xcols(x) x0(0) s(4)
matrix MR_B=e(b)
quietly summarize y
scalar MR_MEAN=r(mean)
capture _rdistnn_messages error, rc(198) message("k() must be positive; got 0.") hint("Set k() to an integer of at least 1.")
test_assert "shared error returns the requested standard code" _rc==198
test_assert "message helper preserves estimation results" e(b)[1,1]==MR_B[1,1]
test_assert "message helper preserves return results" r(mean)==MR_MEAN
capture mata: rdistnn_fail(498,"Saved e(model_groups) is inconsistent. Refit the model.")
test_assert "Mata fatal helper preserves its chosen code" _rc==498
// Plain text logs test host output controls without parsing SMCL formatting.
tempfile message_log
capture program drop MR_outputs
program MR_outputs
    _rdistnn_messages warning, message("ordinary-visible {browse bad} 50%")
    quietly _rdistnn_messages warning, message("quiet-hidden")
    capture _rdistnn_messages error, rc(498) message("capture-hidden")
    capture noisily _rdistnn_messages error, rc(909) message("capture-noisily-visible")
    c_local message_rc=_rc
    capture noisily mata: rdistnn_fail(198,"mata-literal {browse bad} 50%"+char(9)+"control")
    capture noisily MR_quiet_error
end
capture program drop MR_quiet_error
program MR_quiet_error
    quietly _rdistnn_messages error, rc(498) message("quiet-error-visible")
end
quietly log using "`message_log'", text replace name(message_contract)
MR_outputs
quietly log close message_contract
test_assert "capture noisily retains native resource code" `message_rc'==909
file open MR_FILE using "`message_log'", read text
local MR_LOG ""
file read MR_FILE MR_LINE
while r(eof)==0 {
    local MR_LOG `"`MR_LOG' `MR_LINE'"'
    file read MR_FILE MR_LINE
}
file close MR_FILE
mata: MR_LOG=st_local("MR_LOG"); st_numscalar("MR_CHANNELS",strpos(MR_LOG,"ordinary-visible {browse bad} 50%")>0 & strpos(MR_LOG,"quiet-hidden")==0 & strpos(MR_LOG,"capture-hidden")==0 & strpos(MR_LOG,"capture-noisily-visible")>0 & strpos(MR_LOG,"quiet-error-visible")>0)
test_assert "ordinary quiet capture and fatal channels respect Stata controls" MR_CHANNELS
mata: st_numscalar("MR_LITERAL",strpos(MR_LOG,"mata-literal {browse bad} 50% control")>0)
test_assert "Mata fatal messages render dynamic SMCL and percent literally" MR_LITERAL


test_section "bounded dynamic details"
mata {
    MR_LONG=invtokens(J(1,120,"a"),"")
    st_numscalar("MR_BOUNDED",ustrlen(rdistnn_message_value(MR_LONG))==80 & strpos(rdistnn_message_value(MR_LONG),"...")==78)
    st_numscalar("MR_EXAMPLES",rdistnn_message_examples((1::6))=="1, 2, 3, 4, 5 (1 more)" & rdistnn_message_examples((1::5))=="1, 2, 3, 4, 5")
    st_numscalar("MR_CONTROL",rdistnn_message_value("naïve"+char(9)+"{x}%")=="naïve {x}%")
}
test_assert "names are explicitly bounded at 80 characters" MR_BOUNDED
test_assert "five and six locations retain omitted counts" MR_EXAMPLES
test_assert "control characters normalize while Unicode and literal markup survive" MR_CONTROL

test_section "group preparation is reused by public grids and cached consumers"
do "tests/helpers/reporting-instrument.do"
capture program drop _MR_group_reuse_fixture
program define _MR_group_reuse_fixture
    tempfile root
    local MR_SHADOW "`root'_sources"
    mkdir "`MR_SHADOW'"
    capture noisily {
        mata: rdistnn_test_instrument("src/_rdistnn_blocks.ado",st_local("MR_SHADOW")+"/group_count.ado","n=rows(ranks); m=cols(ranks)","st_numscalar("+char(34)+"MR_GROUP_COUNT"+char(34)+",st_numscalar("+char(34)+"MR_GROUP_COUNT"+char(34)+")+1); n=rows(ranks); m=cols(ranks)")
        mata: rdistnn_test_instrument(st_local("MR_SHADOW")+"/group_count.ado",st_local("MR_SHADOW")+"/_rdistnn_blocks.ado","sorted=X[order", "st_numscalar("+char(34)+"MR_DUP_COUNT"+char(34)+",st_numscalar("+char(34)+"MR_DUP_COUNT"+char(34)+")+1); sorted=X[order")
        adopath ++ "`MR_SHADOW'"

        clear
        set obs 24
        gen double x=mod(_n,6)-3
        gen double y=_n^2
        gen double den=1+_n/20
        matrix MR_TARGETS=(0\1\0)
        foreach command in np_dnn cate_dnn clate_dnn {
            local opts ""
            if "`command'"=="cate_dnn" local opts "pseudooutcome"
            if "`command'"=="clate_dnn" local opts "pseudooutcome denominator(den)"
            scalar MR_GROUP_COUNT=0
            scalar MR_DUP_COUNT=0
            quietly `command' y, xcols(x) x0matrix(MR_TARGETS) tuningmode(grid) sgrid(3 4 3 5) sdefault(4) `opts'
            test_assert "`command': one preparation across grid aliases kernels and posting" MR_GROUP_COUNT==1
            test_assert "`command': prepared fitting performs no duplicate-detection sorts" MR_DUP_COUNT==0
            scalar MR_GROUP_COUNT=0
            quietly rdistnn_predict, x0matrix(MR_TARGETS) sefit
            test_assert "`command': identical saved geometry reuses validated groups" MR_GROUP_COUNT==0
            scalar MR_GROUP_COUNT=0
            quietly rdistnn_predict, x0(0.5) sefit
            test_assert "`command': new geometry prepares its own groups once" MR_GROUP_COUNT==1
        }
        gen double mr_d=1
        gen double mr_zero=0
        foreach command in np_tdnn plqr_dnn {
            scalar MR_GROUP_COUNT=0
            scalar MR_DUP_COUNT=0
            if "`command'"=="np_tdnn" quietly np_tdnn y, xcols(x) x0matrix(MR_TARGETS) s1(2) s2(4)
            if "`command'"=="plqr_dnn" quietly plqr_dnn y, x(x) x0(0) d(mr_d) g(mr_zero) md(mr_zero) thetaeval(0) s(4)
            test_assert "`command': weights and covariance reuse one prepared geometry" MR_GROUP_COUNT==1 & MR_DUP_COUNT==0
        }
        scalar MR_GROUP_COUNT=0
        quietly weights_tdnn x, x0matrix(MR_TARGETS) s1(3) s2(5) dim(1)
        test_assert "TDNN weights prepare groups once" MR_GROUP_COUNT==1
        scalar MR_GROUP_COUNT=0
        quietly np_tdnn y, xcols(x) x0matrix(MR_TARGETS) s1(3) s2(5) dim(1)
        test_assert "TDNN estimation reuses weight groups in covariance and posting" MR_GROUP_COUNT==1
    }
    local rc=_rc
    capture adopath - "`MR_SHADOW'"
    quietly do "src/_rdistnn_mata.ado"
    foreach file in group_count _rdistnn_blocks {
        capture erase "`MR_SHADOW'/`file'.ado"
    }
    rmdir "`MR_SHADOW'"
    if `rc' exit `rc'
end
_MR_group_reuse_fixture
program drop _MR_group_reuse_fixture

// Rebind only the host command of the real distance bridge for failures that
// cannot be safely induced in licensed Stata (memory exhaustion and Break).
test_section "native capability absence differs from execution failure"
tempfile MR_DISTANCE_HOST
mata: rdistnn_test_instrument("src/_rdistnn_distance_plugin.ado",st_local("MR_DISTANCE_HOST"),"plugin call rdistnn_plugin, ","_MR_native_host ")
capture program drop _MR_native_host
program _MR_native_host
    version 17
    gettoken operation args : 0
    if inlist("`operation'","distance_api","loo_packed_api") global MR_PROBE_MATRIX `args'
    if "$MR_NATIVE_MODE"=="probe_resource" exit 909
    if "$MR_NATIVE_MODE"=="probe_break" exit 1
    if "$MR_NATIVE_MODE"=="probe_unexpected" exit 498
    if "$MR_NATIVE_MODE"=="absent" exit 198
    if inlist("`operation'","distance_api","loo_packed_api") {
        matrix `args'=(1)
        exit
    }
    gettoken metric args : args
    gettoken X args : args
    gettoken X0 args : args
    gettoken inverse output : args
    global MR_NATIVE_MATRICES `X' `X0' `inverse' `output'
    if "$MR_NATIVE_MODE"=="resource" exit 909
    if "$MR_NATIVE_MODE"=="break" exit 1
    if "$MR_NATIVE_MODE"=="unexpected" exit 498
    // A success code with unchanged missing output is a malformed host result.
    if "$MR_NATIVE_MODE"=="success" matrix `output'=(1\2\3)
end
quietly run "`MR_DISTANCE_HOST'"
global rdistnn_plugin_ok 1
global rdistnn_distance_engine auto
global rdistnn_use_plugin 1
global MR_NATIVE_MODE absent
capture mata: MR_D=rdistnn_distance_matrix((1\2\3),0,"l2")
test_assert "explicit required missing capability fails" _rc==498
macro drop rdistnn_use_plugin
mata: MR_D=rdistnn_distance_plugin((1\2\3),0,"l2",J(0,0,.)); st_numscalar("MR_ABSENT",rows(MR_D)==0)
test_assert "known absent capability returns portable dispatch sentinel" MR_ABSENT
foreach mode in resource probe_resource break unexpected malformed {
    global MR_NATIVE_MODE `mode'
    local expected=498
    if inlist("`mode'","resource","probe_resource") local expected=909
    if "`mode'"=="break" local expected=1
    global rdistnn_use_plugin 1
    capture mata: MR_D=rdistnn_distance_matrix((1\2\3),0,"l2")
    local original_rc=_rc
    test_assert "`mode': actual bridge preserves fatal execution code" `original_rc'==`expected'
    test_assert "`mode': failed native call does not claim Mata execution" "$rdistnn_last_distance_engine"!="mata"
}
// Intercept the real output-copy boundary as well as the host command. The
// production bridge body remains authoritative; no mock bridge is installed.
test_section "native failure skips output retrieval and clears owned matrices"
capture mata: mata drop rdistnn_test_native_copy()
mata:
real matrix rdistnn_test_native_copy(string scalar name, real scalar probe)
{
    string scalar counter
    counter=probe ? "MR_PROBE_COPIES" : "MR_OUTPUT_COPIES"
    st_numscalar(counter,st_numscalar(counter)+1)
    if(st_global("MR_COPY_FAIL")=="1") _error(3900)
    return(st_matrix(name))
}
end
capture program drop _MR_native_copy_fixture
program define _MR_native_copy_fixture
    tempfile host probe output transport_host transport_copy
    do "tests/helpers/e-state.do"
    _rdistnn_test_e_snapshot, name(MR_NATIVE_STATE)
    quietly datasignature
    local caller_data "`r(datasignature)'"
    matrix MR_CALLER_MATRIX=(23,42)
    capture noisily {
        mata: rdistnn_test_instrument("src/_rdistnn_distance_plugin.ado",st_local("host"),"plugin call rdistnn_plugin, ","_MR_native_host ")
        mata: rdistnn_test_instrument(st_local("host"),st_local("probe"),"out=st_matrix(names[1])","out=rdistnn_test_native_copy(names[1],1)")
        // Batch tests leave capability retrieval intact to reach execution.
        mata: rdistnn_test_instrument(st_local("host"),st_local("output"),"out=st_matrix(names[4])","out=rdistnn_test_native_copy(names[4],0)")
        quietly run "`output'"
        foreach mode in success resource break unexpected {
            global MR_NATIVE_MODE `mode'
            foreach copy_fail in 0 1 {
                global MR_COPY_FAIL `copy_fail'
                scalar MR_OUTPUT_COPIES=0
                capture mata: MR_D=rdistnn_distance_plugin((1\2\3),0,"l2_scaled",1)
                local result=_rc
                local expected=498
                if "`mode'"=="resource" local expected=909
                if "`mode'"=="break" local expected=1
                if "`mode'"=="success" local expected=cond(`copy_fail',3900,0)
                test_assert "`mode'/copy`copy_fail': original code survives the output-copy boundary" `result'==`expected'
                test_assert "`mode'/copy`copy_fail': output retrieval runs only on native success" MR_OUTPUT_COPIES==("`mode'"=="success")
                if "`mode'"=="success" & !`copy_fail' {
                    mata: st_numscalar("MR_DISTANCE_VALUES",all(MR_D:==(1\2\3)) & rows(MR_D)==3 & cols(MR_D)==1)
                    test_assert "successful batch keeps exact distances and shape" MR_DISTANCE_VALUES
                }
                // A genuine retrieval failure propagates before this bridge's
                // ordinary cleanup; only native failure/success promise cleanup.
                if "`mode'"!="success" | !`copy_fail' {
                    foreach owned in $MR_PROBE_MATRIX $MR_NATIVE_MATRICES {
                        capture confirm matrix `owned'
                        local exists=_rc==0
                        test_assert "`mode'/copy`copy_fail': owned temporary matrix is cleared" !`exists'
                    }
                }
                foreach owned in $MR_PROBE_MATRIX $MR_NATIVE_MATRICES {
                    capture matrix drop `owned'
                }
                test_assert "native bridge preserves caller matrix" MR_CALLER_MATRIX[1,1]==23 & MR_CALLER_MATRIX[1,2]==42
                quietly datasignature
                local unchanged "`r(datasignature)'"
                test_assert "native bridge preserves caller data" "`unchanged'"=="`caller_data'"
                _rdistnn_test_e_compare, snapshot(MR_NATIVE_STATE)
                local same_e=r(equal)
                test_assert "native bridge preserves complete caller e()" `same_e'
            }
        }
        mata: rdistnn_test_instrument(st_local("probe"),st_local("output"),"out=st_matrix(names[4])","out=rdistnn_test_native_copy(names[4],0)")
        quietly run "`output'"
        foreach mode in absent probe_resource probe_break probe_unexpected {
            global MR_NATIVE_MODE `mode'
            global MR_COPY_FAIL 1
            scalar MR_PROBE_COPIES=0
            scalar MR_OUTPUT_COPIES=0
            capture mata: MR_D=rdistnn_distance_plugin((1\2\3),0,"l2",J(0,0,.))
            local result=_rc
            local expected=498
            if "`mode'"=="absent" local expected=0
            if "`mode'"=="probe_resource" local expected=909
            if "`mode'"=="probe_break" local expected=1
            test_assert "`mode': capability failure/absence skips both copies" MR_PROBE_COPIES==0 & MR_OUTPUT_COPIES==0
            test_assert "`mode': capability preserves its original protocol code" `result'==`expected'
        }
        // The packed-deletion capability bridge has the same host boundary.
        // Arm a failing output read while exercising its actual source body.
        mata: rdistnn_test_instrument("src/_rdistnn_native_transport.ado",st_local("transport_host"),"plugin call rdistnn_plugin, ","_MR_native_host ")
        mata: rdistnn_test_instrument(st_local("transport_host"),st_local("transport_copy"),"result=st_matrix(name)","result=rdistnn_test_native_copy(name,1)")
        quietly run "`transport_copy'"
        foreach mode in absent probe_resource probe_break probe_unexpected success {
            global MR_NATIVE_MODE `mode'
            foreach copy_fail in 0 1 {
                global MR_COPY_FAIL `copy_fail'
                scalar MR_PROBE_COPIES=0
                capture mata: st_numscalar("MR_PACKED_CAPABILITY",rdistnn_native_packed_api())
                local result=_rc
                local expected=498
                if "`mode'"=="absent" local expected=0
                if "`mode'"=="probe_resource" local expected=909
                if "`mode'"=="probe_break" local expected=1
                if "`mode'"=="success" local expected=cond(`copy_fail',3900,0)
                test_assert "packed `mode'/copy`copy_fail': original capability code survives output retrieval" `result'==`expected'
                test_assert "packed `mode'/copy`copy_fail': retrieval runs only after native success" MR_PROBE_COPIES==("`mode'"=="success")
                if !`result' {
                    test_assert "packed `mode': capability value preserves absence/success protocol" MR_PACKED_CAPABILITY==("`mode'"=="success")
                }
                if "`mode'"!="success" | !`copy_fail' {
                    capture confirm matrix $MR_PROBE_MATRIX
                    local exists=_rc==0
                    test_assert "packed `mode'/copy`copy_fail': owned capability matrix is cleared" !`exists'
                }
                capture matrix drop $MR_PROBE_MATRIX
                test_assert "packed capability preserves caller matrix" MR_CALLER_MATRIX[1,1]==23 & MR_CALLER_MATRIX[1,2]==42
                quietly datasignature
                local unchanged "`r(datasignature)'"
                test_assert "packed capability preserves caller data" "`unchanged'"=="`caller_data'"
                _rdistnn_test_e_compare, snapshot(MR_NATIVE_STATE)
                local same_e=r(equal)
                test_assert "packed capability preserves complete caller e()" `same_e'
            }
        }
    }
    local fixture_rc=_rc
    foreach owned in $MR_PROBE_MATRIX $MR_NATIVE_MATRICES {
        capture matrix drop `owned'
    }
    _rdistnn_test_e_snapshot_drop, name(MR_NATIVE_STATE)
    quietly do "src/_rdistnn_mata.ado"
    if `fixture_rc' exit `fixture_rc'
end
_MR_native_copy_fixture
program drop _MR_native_copy_fixture
capture mata: mata drop rdistnn_test_native_copy()
macro drop MR_NATIVE_MODE MR_COPY_FAIL MR_PROBE_MATRIX MR_NATIVE_MATRICES
macro drop rdistnn_use_plugin
macro drop rdistnn_distance_engine
program drop _MR_native_host
quietly do "src/_rdistnn_mata.ado"

test_section "CLATE statuses have one bounded summary per issue family"
clear
set obs 24
gen double x=_n
gen double num=2
gen double den=0
matrix MR_MANY=J(100,1,0)
quietly clate_dnn num, xcols(x) x0matrix(MR_MANY) denominator(den) pseudooutcome s(4) nose
matrix MR_STATUS=e(status)
test_assert "all-unavailable ratios retain all component and status rows" rowsof(MR_STATUS)==100 & rowsof(e(components))==100 & colsof(e(estimates))==100
capture program drop MR_ratio_report
program MR_ratio_report
    _rdistnn_report
end
tempfile ratio_log
quietly log using "`ratio_log'", text replace name(ratio_messages)
MR_ratio_report
quietly log close ratio_messages
file open MR_RATIO_FILE using "`ratio_log'", read text
local MR_UNAVAILABLE=0
local MR_WALD=0
local MR_BOUNDED=0
file read MR_RATIO_FILE MR_LINE
while r(eof)==0 {
    if strpos(`"`MR_LINE'"',"ratio is unavailable") | strpos(`"`MR_LINE'"',"ratios unavailable") local ++MR_UNAVAILABLE
    if strpos(`"`MR_LINE'"',"no pointwise Wald") | strpos(`"`MR_LINE'"',"Wald intervals unavailable") local ++MR_WALD
    if strpos(`"`MR_LINE'"',"95 more") local MR_BOUNDED=1
    file read MR_RATIO_FILE MR_LINE
}
file close MR_RATIO_FILE
test_assert "100 unavailable ratios produce one bounded summary" `MR_UNAVAILABLE'==1 & `MR_BOUNDED'
test_assert "nose does not warn about deliberately absent Wald inference" `MR_WALD'==0

test_section "Mata tie notices consume groups and obey quiet execution"
mata: MR_GROUPS=(1\1\2\3)
capture mata: mata drop rdistnn_duplicate_rows()
mata:
real scalar rdistnn_duplicate_rows(real matrix unused)
{
    rdistnn_fail(498,"Reporting unexpectedly rebuilt equality facts.")
    return(0)
}
end
capture quietly mata: rdistnn_distance_warning(MR_GROUPS)
test_assert "prepared notice calls no duplicate detector or sorting path" _rc==0
capture program drop MR_tie_notices
program MR_tie_notices
    mata: rdistnn_distance_warning(MR_GROUPS)
    quietly mata: rdistnn_distance_warning(MR_GROUPS)
    mata: rdistnn_distance_warning(MR_GROUPS)
end
tempfile tie_log
quietly log using "`tie_log'", text replace name(tie_messages)
MR_tie_notices
quietly log close tie_messages
file open MR_TIE_FILE using "`tie_log'", read text
local MR_NOTICES=0
file read MR_TIE_FILE MR_LINE
while r(eof)==0 {
    if strpos(`"`MR_LINE'"',"warning: Exact distance ties") local ++MR_NOTICES
    file read MR_TIE_FILE MR_LINE
}
file close MR_TIE_FILE
test_assert "two independent ordinary calls warn and the quiet call is silent" `MR_NOTICES'==2
quietly do "src/_rdistnn_mata.ado"

test_section "reordered saved target profiles reuse cached geometry"
clear
set obs 24
gen double x=mod(_n,7)
gen double y=_n
matrix MR_ORDERED=(0\1\2)
quietly np_dnn y, xcols(x) x0matrix(MR_ORDERED) s(4)
matrix MR_REORDERED=(2\0\2\1)
// A deliberately unavailable distance engine distinguishes cache selection
// from a coincidentally equal recomputed result.
global rdistnn_distance_engine invalid
capture quietly rdistnn_predict, x0matrix(MR_REORDERED) sefit
local MR_REORDER_RC=_rc
macro drop rdistnn_distance_engine
test_assert "reordered and repeated retained targets need no distance evaluation" `MR_REORDER_RC'==0

test_section "saved plotting and result notices respect quietly"
quietly cate_dnn y, xcols(x) x0(0) pseudooutcome s(4)
capture program drop MR_quiet_consumers
program MR_quiet_consumers
    quietly rdistnn_results, type(weights)
    quietly cate_dnn_plot, type(weights) name(mr_quiet) replace nodraw
end
tempfile quiet_consumers_log
quietly log using "`quiet_consumers_log'", text replace name(quiet_consumers)
MR_quiet_consumers
quietly log close quiet_consumers
file open MR_Q_FILE using "`quiet_consumers_log'", read text
local MR_Q_NOTICES=0
file read MR_Q_FILE MR_LINE
while r(eof)==0 {
    if strpos(`"`MR_LINE'"',"warning: Exact distance ties") local ++MR_Q_NOTICES
    file read MR_Q_FILE MR_LINE
}
file close MR_Q_FILE
test_assert "quiet saved results and weight plots suppress ordinary tie advisories" `MR_Q_NOTICES'==0

// Capability checks already selected this required weights panel. Its resource
// failure must not turn into a dashboard that silently drops the panel.
capture program drop _cate_dnn_plot_weights
program _cate_dnn_plot_weights
    exit 909
end
capture quietly cate_dnn_plot, type(dashboard) name(mr_fault) replace nodraw
local MR_PANEL_RC=_rc
test_assert "dashboard propagates an unexpected native/resource panel failure" `MR_PANEL_RC'==909
capture program drop _cate_dnn_plot_weights
// Reload actual helpers for later suites without redefining existing programs.
program drop _cate_dnn_plot_estimates _cate_dnn_plot_diagnostics _cate_dnn_plot_tuning _cate_dnn_plot_overlap _cate_dnn_plot_influence _cate_dnn_plot_dashboard _np_dnn_plot_outcome _np_dnn_plot_dashboard _rdistnn_plot_wald_bounds
quietly run "src/_rdistnn_plot_helpers.ado"

test_section "packed native probe rejects malformed capability output"
do "tests/helpers/reporting-instrument.do"
tempfile MR_PACKED_HOST
mata: rdistnn_test_instrument("src/_rdistnn_native_transport.ado",st_local("MR_PACKED_HOST"),"plugin call rdistnn_plugin, ","_MR_packed_host ")
capture program drop _MR_packed_host
program _MR_packed_host
    gettoken operation output : 0
    matrix `output'=(1\1)
end
quietly run "`MR_PACKED_HOST'"
capture mata: MR_PACKED_API=rdistnn_native_packed_api()
local MR_PACKED_RC=_rc
test_assert "successful malformed packed capability shape is fatal" `MR_PACKED_RC'==498
program drop _MR_packed_host
quietly do "src/_rdistnn_mata.ado"

test_section "shared public notices cover grids weak localization replay and reuse"
do "tests/helpers/reporting-instrument.do"
clear
set obs 24
gen double x=mod(_n,6)
gen double y=_n
gen double den=0
capture program drop MR_grid_notice
program MR_grid_notice
    np_dnn y, xcols(x) x0(0) tuningmode(grid) sgrid(3 30) sdefault(3)
end
tempfile grid_log
quietly log using "`grid_log'", text replace name(grid_notice)
MR_grid_notice
quietly log close grid_notice
_rdistnn_test_notice_count using "`grid_log'", pattern("infeasible optional")
test_assert "infeasible optional candidates have one bounded warning and retain their row" r(count)==1 & rowsof(e(tuning_candidates))==2 & e(tuning_candidates)[2,4]==0

quietly np_dnn y, xcols(x) x0(0) s(4) k(4)
capture program drop MR_replay_notice
program MR_replay_notice
    np_dnn
end
tempfile replay_log
quietly log using "`replay_log'", text replace name(replay_notice)
MR_replay_notice
quietly log close replay_notice
_rdistnn_test_notice_count using "`replay_log'", pattern("warning: Exact distance ties")
test_assert "an independent estimator replay emits its tie advisory once" r(count)==1
_rdistnn_test_notice_count using "`replay_log'", pattern("weak localization")
test_assert "uniform DNN weights warn at the existing .9*n ESS guide" r(count)==1

quietly clate_dnn y, xcols(x) x0(0) denominator(den) pseudooutcome s(4)
capture program drop MR_reuse_notice
program MR_reuse_notice
    rdistnn_predict, x0(1) sefit
end
tempfile reuse_log
quietly log using "`reuse_log'", text replace name(reuse_notice)
MR_reuse_notice
quietly log close reuse_notice
_rdistnn_test_notice_count using "`reuse_log'", pattern("ratios unavailable")
test_assert "ordinary prediction summarizes unavailable ratios once" r(count)==1
_rdistnn_test_notice_count using "`reuse_log'", pattern("Wald intervals unavailable")
test_assert "requested prediction uncertainty summarizes unavailable intervals once" r(count)==1

quietly replace x=0
quietly np_tdnn y, xcols(x) x0(0) s1(1) s2(2) dim(1)
capture program drop MR_tdnn_notice
program MR_tdnn_notice
    np_tdnn
end
tempfile tdnn_log
quietly log using "`tdnn_log'", text replace name(tdnn_notice)
MR_tdnn_notice
local MR_TDNN_WEAK_COUNT=r(n_weak_targets)
quietly log close tdnn_notice
test_assert "TDNN replay returns its existing near-uniform signed concentration count" e(n_weak_targets)==1 & `MR_TDNN_WEAK_COUNT'==1
_rdistnn_test_notice_count using "`tdnn_log'", pattern("signed-weight concentration")
test_assert "TDNN near-uniform concentration uses signed-weight terminology" r(count)==1

quietly cate_dnn y, xcols(x) x0(0) pseudooutcome s(4)
capture program drop MR_quiet_extra_views
program MR_quiet_extra_views
    quietly rdistnn_results, type(influence)
    quietly cate_dnn_plot, type(weights) view(kernel) name(mr_quiet_kernel) replace nodraw
    quietly cate_dnn_plot, type(influence) name(mr_quiet_influence) replace nodraw
end
tempfile quiet_views_log
quietly log using "`quiet_views_log'", text replace name(quiet_extra_views)
MR_quiet_extra_views
quietly log close quiet_extra_views
_rdistnn_test_notice_count using "`quiet_views_log'", pattern("warning: Exact distance ties")
test_assert "quiet influence and kernel consumers suppress text advisories" r(count)==0

test_section "ordinary plot calls emit each family once across panels"
quietly np_dnn y, xcols(x) x0(0) s(4)
capture program drop MR_ordinary_plot
program MR_ordinary_plot
    args view
    np_dnn_plot, type(`view') name(mr_ordinary_`view') replace nodraw
end
foreach view in weights dashboard {
    tempfile plot_notice_log
    quietly log using "`plot_notice_log'", text replace name(ordinary_plot)
    MR_ordinary_plot `view'
    quietly log close ordinary_plot
    _rdistnn_test_notice_count using "`plot_notice_log'", pattern("warning: Exact distance ties")
    test_assert "ordinary `view' plot emits one tie advisory" r(count)==1
}

test_section "notice concentration comes from validated diagnostic rows"
quietly np_dnn y, xcols(x) x0(0) k(4) s(4)
capture program drop MR_bad_concentration
program MR_bad_concentration, eclass
    tempname bad
    matrix `bad'=e(diagnostics)
    matrix `bad'[1,1]=.
    ereturn matrix diagnostics=`bad', copy
end
MR_bad_concentration
capture quietly np_dnn
local bad_rc=_rc
test_assert "corrupt cached concentration fails before notice consumption" `bad_rc'==498
program drop MR_bad_concentration

test_section "singleton nose fits keep the shared localization guard"
clear
set obs 1
gen double x=0
gen double y=2
capture program drop MR_singleton
program MR_singleton
    args engine
    global rdistnn_distance_engine `engine'
    np_dnn y, xcols(x) x0(0) s(1) nose
end
foreach engine in mata plugin {
    tempfile singleton_log
    quietly log using "`singleton_log'", text replace name(singleton_notice)
    MR_singleton `engine'
    local MR_SINGLETON_TIE_COUNT=r(n_tied_targets)
    local MR_SINGLETON_WEAK_COUNT=r(n_weak_targets)
    quietly log close singleton_notice
    test_assert "`engine': singleton notice flags and counts stay zero" e(n_tied_targets)==0 & e(n_weak_targets)==0 & `MR_SINGLETON_TIE_COUNT'==0 & `MR_SINGLETON_WEAK_COUNT'==0
    _rdistnn_test_notice_count using "`singleton_log'", pattern("weak localization")
    test_assert "`engine': singleton nose skips near-uniform localization warning" r(count)==0
}
macro drop rdistnn_distance_engine


test_section "quiet notice facts are machine-readable"
clear
set obs 12
gen double x=mod(_n,3)
gen double y=_n
matrix MR_FACT_TARGETS=(0\1\0)
quietly np_dnn y, xcols(x) x0matrix(MR_FACT_TARGETS) s(1) nose
local MR_FIELDS : e(matrices)
local MR_HAS_FACTS : list posof "target_notices" in MR_FIELDS
test_assert "fresh quiet fits retain per-target notice facts" `MR_HAS_FACTS'>0
test_assert "fresh quiet fits retain complete notice counts" e(n_tied_targets)==3 & e(n_weak_targets)==3
test_assert "fresh fit returns notice counts" r(n_tied_targets)==3 & r(n_weak_targets)==3
quietly rdistnn_results, targets(3 1 3)
test_assert "selected result counts deduplicate target IDs" r(n_tied_targets)==2 & r(n_weak_targets)==2
quietly rdistnn_predict, x0matrix(MR_FACT_TARGETS)
test_assert "derived predictions return their batch counts" r(n_tied_targets)==3 & r(n_weak_targets)==3
quietly weights_dnn x, x0matrix(MR_FACT_TARGETS) s(1)
test_assert "standalone weights return tie and concentration counts" r(n_tied_targets)==3 & r(n_weak_targets)==3


test_section "derived notice failures restore the original estimation state"
do "tests/helpers/e-state.do"
quietly np_dnn y, xcols(x) x0(0) s(2)
_rdistnn_test_e_snapshot, name(MR_NOTICE_STATE)
findfile "_rdistnn_report_notices.ado"
local MR_NOTICE_SOURCE "`r(fn)'"
program drop _rdistnn_report_notices
program define _rdistnn_report_notices, rclass
    syntax [, TARGETS(numlist integer >=1)]
    if c(noisily) exit $MR_NOTICE_RC
    return scalar n_tied_targets=e(n_tied_targets)
    return scalar n_weak_targets=e(n_weak_targets)
end
foreach fail_rc in 1 909 498 {
    global MR_NOTICE_RC `fail_rc'
    capture noisily rdistnn_predict, x0(0.5)
    local MR_NOTICE_RESULT=_rc
    test_assert "notice failure `fail_rc' preserves its return code" `MR_NOTICE_RESULT'==`fail_rc'
    _rdistnn_test_e_compare, snapshot(MR_NOTICE_STATE)
    test_assert "notice failure `fail_rc' restores complete original e()" r(equal)
}
macro drop MR_NOTICE_RC
program drop _rdistnn_report_notices
quietly run "`MR_NOTICE_SOURCE'"
_rdistnn_test_e_snapshot_drop, name(MR_NOTICE_STATE)

test_section "fatal geometry messages stay clean and literal"
capture program drop MR_geometry_error
program define MR_geometry_error
    capture noisily mata: MR_BAD=rdistnn_fit_geometry((1\2),"l2_scaled","unknown",J(1,0,.),1e-12)
    c_local MR_GEOMETRY_RC=_rc
end
tempfile geometry_log
quietly log using "`geometry_log'", text replace name(geometry_messages)
MR_geometry_error
quietly log close geometry_messages
file open MR_G_FILE using "`geometry_log'", read text
local MR_GEOMETRY_TRACE=0
local MR_GEOMETRY_MESSAGE=0
file read MR_G_FILE MR_LINE
while r(eof)==0 {
    if strpos(`"`MR_LINE'"',"<istmt>") | strpos(`"`MR_LINE'"',"r(198)") | strpos(`"`MR_LINE'"',"invalid syntax") local MR_GEOMETRY_TRACE=1
    if strpos(`"`MR_LINE'"',"scalerule() must be sd, mad, iqr, range, or custom.") local MR_GEOMETRY_MESSAGE=1
    file read MR_G_FILE MR_LINE
}
file close MR_G_FILE
test_assert "geometry option error preserves standard syntax code" `MR_GEOMETRY_RC'==198
test_assert "geometry fatal has problem text without duplicate generic error or traceback" `MR_GEOMETRY_MESSAGE' & !`MR_GEOMETRY_TRACE'
program drop MR_geometry_error

test_section "all-target notice counts avoid numlist expansion limits"
_rdistnn_test_e_snapshot, name(MR_LARGE_NOTICE_STATE)
tempname MR_LARGE_ORIGINAL
_estimates hold `MR_LARGE_ORIGINAL', copy restore
capture program drop MR_large_notice_state
program define MR_large_notice_state, eclass
    tempname facts
    matrix `facts'=J(3000,2,0)
    matrix `facts'[3000,1]=1
    matrix `facts'[1,2]=1
    ereturn clear
    ereturn matrix target_notices=`facts', copy
    ereturn scalar m=3000
    ereturn scalar N=12
    ereturn local cmd "np_dnn"
end
MR_large_notice_state
capture quietly _rdistnn_report_notices
local MR_LARGE_RC=_rc
local MR_LARGE_TIES=.
local MR_LARGE_WEAK=.
if !`MR_LARGE_RC' {
    local MR_LARGE_TIES=r(n_tied_targets)
    local MR_LARGE_WEAK=r(n_weak_targets)
}
_estimates unhold `MR_LARGE_ORIGINAL'
test_assert "all-target notices accept batches above the numlist limit" `MR_LARGE_RC'==0
test_assert "large-batch notices count the whole retained matrix" `MR_LARGE_TIES'==1 & `MR_LARGE_WEAK'==1
_rdistnn_test_e_compare, snapshot(MR_LARGE_NOTICE_STATE)
test_assert "large-batch notice fixture restores complete caller e()" r(equal)
_rdistnn_test_e_snapshot_drop, name(MR_LARGE_NOTICE_STATE)
program drop MR_large_notice_state

test_section "notice log counter treats captured output as literal data"
capture program drop MR_literal_logged_lines
program define MR_literal_logged_lines
    // Log wrapping can separate opening and closing compound quotes.
    mata: display("notice-literal "+char(34)+"double"+char(34)+" "+char(96)+char(34)+"compound {browse bad} 50%",1)
    mata: display("notice-literal compound"+char(34)+char(39)+" {browse bad} 50%",1)
end
tempfile MR_LITERAL_NOTICE_LOG
quietly log using "`MR_LITERAL_NOTICE_LOG'", text replace name(literal_notice_log)
MR_literal_logged_lines
quietly log close literal_notice_log
capture noisily _rdistnn_test_notice_count using "`MR_LITERAL_NOTICE_LOG'", pattern("notice-literal")
local MR_LITERAL_COUNT_RC=_rc
local MR_LITERAL_COUNT=.
if !`MR_LITERAL_COUNT_RC' local MR_LITERAL_COUNT=r(count)
test_assert "log counter accepts double compound quotes SMCL and percent" `MR_LITERAL_COUNT_RC'==0
test_assert "literal log counter retains the exact occurrence count" `MR_LITERAL_COUNT'==2
program drop MR_literal_logged_lines
