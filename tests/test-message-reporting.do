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
local MR_ROOT : pwd
local MR_ROOT = subinstr("`MR_ROOT'", "\", "/", .)
local MR_SHADOW "`MR_ROOT'/build/reporting-implementation-20261022/instrument"
capture mkdir "build"
capture mkdir "build/reporting-implementation-20261022"
capture mkdir "`MR_SHADOW'"
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
adopath - "`MR_SHADOW'"
quietly do "src/_rdistnn_mata.ado"

// Rebind only the host command of the real distance bridge for failures that
// cannot be safely induced in licensed Stata (memory exhaustion and Break).
test_section "native capability absence differs from execution failure"
mata: rdistnn_test_instrument("src/_rdistnn_distance_plugin.ado",st_local("MR_SHADOW")+"/distance-host.ado","plugin call rdistnn_plugin, ","_MR_native_host ")
capture program drop _MR_native_host
program _MR_native_host
    version 17
    gettoken operation args : 0
    if "$MR_NATIVE_MODE"=="probe_resource" exit 909
    if "$MR_NATIVE_MODE"=="absent" exit 198
    if "`operation'"=="distance_api" {
        matrix `args'=(1)
        exit
    }
    if "$MR_NATIVE_MODE"=="resource" exit 909
    if "$MR_NATIVE_MODE"=="break" exit 1
    if "$MR_NATIVE_MODE"=="unexpected" exit 498
    // A success code with unchanged missing output is a malformed host result.
end
quietly run "`MR_SHADOW'/distance-host.ado"
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
macro drop MR_NATIVE_MODE
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
local MR_ROOT : pwd
local MR_ROOT=subinstr("`MR_ROOT'","\","/",.)
local MR_SHADOW "`MR_ROOT'/build/reporting-implementation-20261022/instrument"
mata: rdistnn_test_instrument("src/_rdistnn_native_transport.ado",st_local("MR_SHADOW")+"/packed-host.ado","plugin call rdistnn_plugin, ","_MR_packed_host ")
capture program drop _MR_packed_host
program _MR_packed_host
    gettoken operation output : 0
    matrix `output'=(1\1)
end
quietly run "`MR_SHADOW'/packed-host.ado"
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
quietly log close tdnn_notice
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
    quietly log close singleton_notice
    _rdistnn_test_notice_count using "`singleton_log'", pattern("weak localization")
    test_assert "`engine': singleton nose skips near-uniform localization warning" r(count)==0
}
macro drop rdistnn_distance_engine
