// Optional native loading and saved geometry distinguish absence from failure.
version 17
do "tests/test_framework.do"
do "tests/helpers/reporting-instrument.do"
do "tests/helpers/e-state.do"
do "tests/helpers/caller-data.do"
capture mkdir "build"
capture mkdir "build/test-optional-boundaries"
local shadow "build/test-optional-boundaries"
local old_distance "$rdistnn_distance_engine"
local old_jackknife "$rdistnn_jackknife_engine"
local old_use_plugin "$rdistnn_use_plugin"
global rdistnn_distance_engine mata
global rdistnn_jackknife_engine mata
global rdistnn_use_plugin 0

test_section "optional plugin declaration preserves unexpected return codes"
capture program drop _OB_plugin_declare
program define _OB_plugin_declare
    args rc
    exit `rc'
end
foreach fail_rc in 1 909 498 198 {
    mata: rdistnn_test_instrument("src/_rdistnn_mata.ado",st_local("shadow")+"/plugin-failure.ado","program define rdistnn_plugin, plugin using(rdistnn.plugin)","_OB_plugin_declare "+st_local("fail_rc"))
    macro drop rdistnn_plugin_ok
    capture noisily quietly run "`shadow'/plugin-failure.ado"
    local result_rc=_rc
    test_assert "plugin declaration failure `fail_rc' retains its return code" `result_rc'==`fail_rc'
    test_assert "plugin declaration failure `fail_rc' does not cache unavailability" "$rdistnn_plugin_ok"==""
    test_assert "plugin declaration failure `fail_rc' preserves engine selectors" "$rdistnn_distance_engine"=="mata" & "$rdistnn_jackknife_engine"=="mata" & "$rdistnn_use_plugin"=="0"
}
program drop _OB_plugin_declare

test_section "actual missing, nonloadable and successful plugin declarations"
capture program drop rdistnn_plugin
mata: rdistnn_test_instrument("src/_rdistnn_mata.ado",st_local("shadow")+"/plugin-missing.ado","using(rdistnn.plugin)","using(rdistnn_absent_boundary_probe.plugin)")
macro drop rdistnn_plugin_ok
capture noisily quietly run "`shadow'/plugin-missing.ado"
local missing_rc=_rc
test_assert "an absent optional plugin selects portable capability" `missing_rc'==0 & "$rdistnn_plugin_ok"=="0"

tempfile invalid_plugin
tempname handle
file open `handle' using "`invalid_plugin'", write text replace
file write `handle' "not a native shared library" _n
file close `handle'
mata: rdistnn_test_instrument("src/_rdistnn_mata.ado",st_local("shadow")+"/plugin-unloadable.ado","using(rdistnn.plugin)","using("+char(34)+st_local("invalid_plugin")+char(34)+")")
macro drop rdistnn_plugin_ok
capture noisily quietly run "`shadow'/plugin-unloadable.ado"
local unloadable_rc=_rc
test_assert "a nonloadable optional binary selects portable capability" `unloadable_rc'==0 & "$rdistnn_plugin_ok"=="0"

// Native availability is determined by the actual host binary, not a test flag.
capture program define OB_NATIVE_CONTROL, plugin using(rdistnn.plugin)
local declaration_rc=_rc
local native_available=inlist(`declaration_rc',0,110)
test_assert "host declaration has a supported success or unavailable outcome" inlist(`declaration_rc',0,110,601,9999)
capture program drop OB_NATIVE_CONTROL
macro drop rdistnn_plugin_ok
quietly run "src/_rdistnn_mata.ado"
test_assert "ordinary loading matches the actual host native capability" real("$rdistnn_plugin_ok")==`native_available'
macro drop rdistnn_plugin_ok
quietly run "src/_rdistnn_mata.ado"
test_assert "already-declared native handles remain available on reload" real("$rdistnn_plugin_ok")==`native_available'

test_section "TDNN optional geometry uses returned matrix membership"
clear
set obs 12
gen double ob_x1=_n
gen double ob_x2=mod(_n,3)
gen double ob_y=1+.5*ob_x1-.2*ob_x2
matrix OB_TARGETS=(.2,.1 \ 8.7,.3)
quietly np_tdnn ob_y, xcols(ob_x1 ob_x2) x0matrix(OB_TARGETS) s1(3) s2(6) nose
matrix OB_EXPECTED=e(estimates)
// Omit only the optional return from the real weights implementation.
mata: rdistnn_test_instrument("src/weights_tdnn.ado",st_local("shadow")+"/weights-without-geometry.ado","if "+char(96)+"geometry_rc"+char(39)+"==0 return matrix geometry=","if 0 return matrix geometry=")
program drop weights_tdnn
quietly run "`shadow'/weights-without-geometry.ado"
capture noisily quietly np_tdnn ob_y, xcols(ob_x1 ob_x2) x0matrix(OB_TARGETS) s1(3) s2(6) nose
local absent_geometry_rc=_rc
local fields : e(matrices)
local geometry_present : list posof "model_geometry" in fields
test_assert "unscaled TDNN permits actual absence of optional geometry" `absent_geometry_rc'==0
test_assert "absent geometry does not become a fabricated missing matrix" `geometry_present'==0
test_assert "omitting optional unscaled geometry preserves estimates" mreldif(e(estimates),OB_EXPECTED)==0
program drop weights_tdnn
quietly run "src/weights_tdnn.ado"

quietly weights_tdnn ob_x1 ob_x2, x0matrix(OB_TARGETS) s1(3) s2(6) distance(l2_scaled)
matrix OB_GEOMETRY=r(geometry)
quietly np_tdnn ob_y, xcols(ob_x1 ob_x2) x0matrix(OB_TARGETS) s1(3) s2(6) distance(l2_scaled) nose
test_assert "returned frozen geometry is copied exactly on success" mreldif(e(model_geometry),OB_GEOMETRY)==0

test_section "TDNN geometry copy failures propagate through caller cleanup"
capture program drop _OB_geometry_copy
program define _OB_geometry_copy
    args target rc
    // A partially produced temporary destination must also be cleaned up.
    matrix `target'=(987)
    exit `rc'
end
tempfile caller_data
quietly save "`caller_data'", replace
_rdistnn_test_data_snapshot
foreach fail_rc in 1 909 498 198 601 {
    quietly np_tdnn ob_y, xcols(ob_x1 ob_x2) x0matrix(OB_TARGETS) s1(3) s2(6) distance(l2_scaled) nose
    _rdistnn_test_e_snapshot, name(OB_STATE)
    local matrices_before : all matrices
    mata: rdistnn_test_instrument("src/np_tdnn.ado",st_local("shadow")+"/geometry-failure.ado","matrix "+char(96)+"geom"+char(39)+"=r(geometry)","_OB_geometry_copy "+char(96)+"geom"+char(39)+" "+st_local("fail_rc"))
    capture program drop np_tdnn
    capture program drop _np_tdnn_Display
    quietly run "`shadow'/geometry-failure.ado"
    capture noisily quietly np_tdnn ob_y, xcols(ob_x1 ob_x2) x0matrix(OB_TARGETS) s1(3) s2(6) distance(l2_scaled) nose
    local copy_rc=_rc
    test_assert "geometry copy failure `fail_rc' retains its return code" `copy_rc'==`fail_rc'
    _rdistnn_test_e_compare, snapshot(OB_STATE)
    test_assert "geometry copy failure `fail_rc' preserves complete prior estimates" r(equal)
    _rdistnn_test_data_compare, data("`caller_data'")
    test_assert "geometry copy failure `fail_rc' preserves caller data" r(equal)
    local matrices_after : all matrices
    local extra : list matrices_after - matrices_before
    local absent : list matrices_before - matrices_after
    test_assert "geometry copy failure `fail_rc' cleans temporary matrices" "`extra'`absent'"==""
    _rdistnn_test_e_snapshot_drop, name(OB_STATE)
    capture program drop np_tdnn
    capture program drop _np_tdnn_Display
    quietly run "src/np_tdnn.ado"
}
program drop _OB_geometry_copy

test_section "optional matrix and program probes preserve nonabsence failures"
capture program drop _OB_probe
program define _OB_probe
    exit $OB_PROBE_RC
end
capture program drop _OB_probe_fixture
program define _OB_probe_fixture
    tempfile injected caller_data graph_file
    quietly save "`caller_data'", replace
    _rdistnn_test_data_snapshot
    capture noisily {
        foreach code in 1 909 498 {
            quietly np_dnn ob_y, xcols(ob_x1 ob_x2) x0matrix(OB_TARGETS) s(3)
            _rdistnn_test_e_snapshot, name(OB_PROBE_STATE)
            global OB_PROBE_RC `code'
            mata: rdistnn_test_instrument("src/weights_tdnn.ado",st_local("injected"),"capture confirm matrix "+char(96)+"geom"+char(39),"capture _OB_probe")
            program drop weights_tdnn
            quietly run "`injected'"
            capture quietly weights_tdnn ob_x1 ob_x2, x0(.2 .1) s1(3) s2(6) generate(ob_probe_output)
            local result=_rc
            test_assert "optional matrix probe preserves native`code'" `result'==`code'
            capture confirm variable ob_probe_output
            local output_rc=_rc
            test_assert "failed geometry probe creates no output variable" `output_rc'==111
            program drop weights_tdnn
            quietly run "src/weights_tdnn.ado"

            foreach family in np cate clate {
                mata: rdistnn_test_instrument("src/"+st_local("family")+"_dnn_plot.ado",st_local("injected"),"capture program list _rdistnn_plot_dispatch","capture _OB_probe")
                capture program drop `family'_dnn_plot
                quietly run "`injected'"
                capture quietly `family'_dnn_plot, nodraw
                local result=_rc
                test_assert "`family': dispatcher probe preserves native`code'" `result'==`code'
                program drop `family'_dnn_plot
                quietly run "src/`family'_dnn_plot.ado"
            }
            mata: rdistnn_test_instrument("src/_rdistnn_plot_dispatch.ado",st_local("injected"),"capture program list _cate_dnn_plot_estimates","capture _OB_probe")
            foreach helper in _rdistnn_plot_dispatch _rdistnn_plot_level _rdistnn_plot_candidates _rdistnn_plot_sample_guard _rdistnn_plot_export _rdistnn_plot_add_panel {
                capture program drop `helper'
            }
            quietly run "`injected'"
            capture quietly np_dnn_plot, nodraw
            local result=_rc
            test_assert "plot helper probe preserves native`code'" `result'==`code'
            program drop _rdistnn_plot_dispatch _rdistnn_plot_level _rdistnn_plot_candidates _rdistnn_plot_sample_guard _rdistnn_plot_export _rdistnn_plot_add_panel
            quietly run "src/_rdistnn_plot_dispatch.ado"
            _rdistnn_test_e_compare, snapshot(OB_PROBE_STATE)
            test_assert "optional probe failures preserve complete estimates" r(equal)
            _rdistnn_test_data_compare, data("`caller_data'")
            test_assert "optional probe failures preserve caller data" r(equal)
            _rdistnn_test_e_snapshot_drop, name(OB_PROBE_STATE)

            quietly np_tdnn ob_y, xcols(ob_x1 ob_x2) x0matrix(OB_TARGETS) s1(3) s2(6)
            _rdistnn_test_e_snapshot, name(OB_PROBE_STATE)
            mata: rdistnn_test_instrument("src/np_tdnn_plot.ado",st_local("injected"),"capture program list ","capture _OB_probe ")
            capture program drop np_tdnn_plot
            quietly run "`injected'"
            foreach view in estimates weights {
                local saving ""
                if "`view'"=="weights" local saving `"saving(`graph_file'.gph)"'
                capture quietly np_tdnn_plot, type(`view') name(ob_probe) replace nodraw `saving'
                local result=_rc
                test_assert "TDNN `view': helper probe preserves native`code'" `result'==`code'
                _rdistnn_test_e_compare, snapshot(OB_PROBE_STATE)
                test_assert "TDNN probe failure preserves complete estimates" r(equal)
                _rdistnn_test_data_compare, data("`caller_data'")
                test_assert "TDNN probe failure preserves caller data" r(equal)
            }
            program drop np_tdnn_plot
            quietly run "src/np_tdnn_plot.ado"
            _rdistnn_test_e_snapshot_drop, name(OB_PROBE_STATE)
        }
    }
    local fixture_rc=_rc
    foreach command in weights_tdnn np_dnn_plot cate_dnn_plot clate_dnn_plot np_tdnn_plot {
        capture program drop `command'
        quietly run "src/`command'.ado"
    }
    foreach helper in _rdistnn_plot_dispatch _rdistnn_plot_level _rdistnn_plot_candidates _rdistnn_plot_sample_guard _rdistnn_plot_export _rdistnn_plot_add_panel {
        capture program drop `helper'
    }
    quietly run "src/_rdistnn_plot_dispatch.ado"
    _rdistnn_test_e_snapshot_drop, name(OB_PROBE_STATE)
    macro drop OB_PROBE_RC
    if `fixture_rc' exit `fixture_rc'
end
_OB_probe_fixture
program drop _OB_probe_fixture _OB_probe
matrix drop OB_TARGETS OB_EXPECTED OB_GEOMETRY
global rdistnn_distance_engine "`old_distance'"
global rdistnn_jackknife_engine "`old_jackknife'"
global rdistnn_use_plugin "`old_use_plugin'"
