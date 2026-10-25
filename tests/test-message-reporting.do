// Public error channels, standard codes, bounded notices and caller state.
// Geometry-work counters and native/plot fault injection have separate owners;
// this suite keeps cold validation and observable notice behavior together.
version 17
do "tests/test_framework.do"
test_section "shared reporting channels and state"
clear
set obs 12
gen double x=_n
gen double y=2*x
// Exercise errors before any Mata core or plot auxiliary program has loaded.
discard
mata: mata clear
do "tests/test_framework.do"
capture weights_dnn x, x0(0) method(unsupported)
local cold_option_rc=_rc
test_assert "cold custom option error retains198" `cold_option_rc'==198
capture np_dnn_plot
local cold_plot_rc=_rc
test_assert "cold plot error retains missing-estimates301" `cold_plot_rc'==301
capture program list MR_nonexistent_helper
local missing_program_rc=_rc
test_assert "native missing program probe is111" `missing_program_rc'==111
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
    mata: st_local("literal", "ado-literal "+char(34)+"quoted"+char(34)+" {browse bad} 50%"+char(9)+"tab"+char(13)+"cr"+char(10)+"lf")
    capture noisily _rdistnn_messages error, rc(198) message(`"`macval(literal)'"')
    c_local literal_rc=_rc
    mata: st_local("compound", "ado-compound "+char(96)+char(34)+"literal {browse bad} 50%"+char(34)+char(39))
    capture noisily _rdistnn_messages error, rc(198) message(`"`macval(compound)'"')
    c_local compound_rc=_rc
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
mata: st_numscalar("MR_ADO_LITERAL",strpos(MR_LOG,"ado-literal "+char(34)+"quoted"+char(34)+" {browse bad} 50% tab cr lf")>0)
test_assert "ado fatal messages preserve quotes and normalize control characters" MR_ADO_LITERAL & `literal_rc'==198
mata: st_numscalar("MR_ADO_COMPOUND",strpos(MR_LOG,"ado-compound "+char(96)+char(34)+"literal {browse bad} 50%"+char(34)+char(39))>0)
test_assert "ado fatal messages preserve balanced compound quote text" MR_ADO_COMPOUND & `compound_rc'==198


test_section "standard native codes and literal method messages"
foreach family in np_dnn cate_dnn clate_dnn {
    local extra ""
    if "`family'"=="cate_dnn" local extra "pseudooutcome"
    if "`family'"=="clate_dnn" local extra "pseudooutcome denominator(y)"
    capture `family' y, xcols(x) x0matrix(MR_matrix_does_not_exist) s(4) `extra'
    local MR_rc=_rc
    test_assert "`family': missing target matrix retains native111" `MR_rc'==111
}
capture mata: MR_bad=rdistnn_weighted_signals(J(2,1,1),J(3,1,1))
test_assert "weighted signal shape mismatch uses standard503" _rc==503
capture mata: MR_bad=rdistnn_jk_cov_from_ranks((1::3),J(2,1,1),1,2,"exact")
test_assert "covariance shape mismatch uses standard503" _rc==503

tempfile MR_literal_log
capture program drop MR_method_literal
program MR_method_literal
    capture noisily mata: MR_bad=rdistnn_weights_from_ranks((1::3),1,2,"{bf:bad} 50%")
end
quietly log using "`MR_literal_log'", text replace name(MR_literal_contract)
quietly MR_method_literal
quietly log close MR_literal_contract
mata: MR_literal_lines=cat(st_local("MR_literal_log")); st_numscalar("MR_LITERAL_METHOD",any(strpos(MR_literal_lines,"unknown method '{bf:bad} 50%'"):>0))
test_assert "actual Mata method error treats dynamic SMCL as literal data" MR_LITERAL_METHOD
program drop MR_method_literal

test_section "native confirmation codes and pure shape guards"
foreach family in np_dnn cate_dnn clate_dnn {
    local extra ""
    if "`family'"=="cate_dnn" local extra "pseudooutcome"
    if "`family'"=="clate_dnn" local extra "pseudooutcome denominator(y)"
    matrix MR_WRONG_WIDTH=(0,0)
    capture `family' y, xcols(x) x0matrix(MR_WRONG_WIDTH) s(4) `extra'
    local MR_rc=_rc
    test_assert "`family': incompatible target columns use503" `MR_rc'==503
}
foreach command in weights_tdnn np_tdnn {
    local vars "x"
    local opts ""
    if "`command'"=="np_tdnn" {
        local vars "y"
        local opts "xcols(x)"
    }
    capture `command' `vars', `opts' x0(0) s1(2) s2(4) dim(1.5)
    local MR_rc=_rc
    test_assert "`command': integer confirmation preserves native7" `MR_rc'==7
}
capture mata: MR_bad=rdistnn_distance_groups((1::3),J(2,1,0))
test_assert "distance group alignment uses503" _rc==503
quietly weights_tdnn x, x0(0) s1(2) s2(4)
capture mata: MR_bad=rdistnn_tdnn_diagnostics(J(3,1,1),J(2,1,0))
test_assert "signed diagnostic alignment uses503" _rc==503
capture mata: MR_bad=rdistnn_weighted_quantile(J(0,1,.),J(0,1,.),.5)
test_assert "empty quantiles use standard no-observations2000" _rc==2000

test_section "bounded dynamic details"
mata {
    MR_LONG=invtokens(J(1,120,"a"),"")
    st_numscalar("MR_BOUNDED",ustrlen(rdistnn_message_value(MR_LONG))==80 & strpos(rdistnn_message_value(MR_LONG),"...")==78)
    st_numscalar("MR_CONTROL",rdistnn_message_value("naïve"+char(9)+"{x}%")=="naïve {x}%")
}
test_assert "names are explicitly bounded at 80 characters" MR_BOUNDED
test_assert "control characters normalize while Unicode and literal markup survive" MR_CONTROL

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

clear
set obs 24
gen double x=mod(_n,7)
gen double y=_n
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
