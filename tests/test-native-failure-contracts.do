// Native capability, execution and output-copy failures are distinct contracts.
// Instrument actual bridge source at host/copy boundaries, preserve original
// codes and caller state, then reload production providers. Injected Break
// proves code propagation; it does not establish live interruption latency.
version 17
do "tests/test_framework.do"
do "tests/helpers/reporting-instrument.do"
clear
set obs 12
gen double x=_n
gen double y=2*x
quietly np_dnn y, xcols(x) x0(0) s(4)
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
