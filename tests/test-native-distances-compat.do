// Run ONLY in a separate process whose rdistnn_plugin declaration points to
// build/native/baseline/rdistnn.plugin (pre-distance-API supported binary).
version 17
// The ordinary suite keeps its current native module loaded; do not drop or
// replace that module. Parent runtime driver defines this global explicitly.
if "$rdistnn_test_old_distance_plugin"!="1" exit 0
quietly do "src/_rdistnn_mata.ado"
do "tests/test_framework.do"
test_section "native distances: real old native capability and failed operation"
// This substantial batch qualifies for automatic native dispatch, so the
// fallback assertions below exercise real API/operation failures.
mata: st_matrix("NC_X",J(5000,16,0));st_matrix("NC_T",J(256,16,1))
global rdistnn_jackknife_engine mata
global rdistnn_plugin_ok 1
capture macro drop rdistnn_use_plugin
global rdistnn_distance_engine mata
mata: nc_ref=rdistnn_distance_matrix(st_matrix("NC_X"),st_matrix("NC_T"),"l2")
// Empty capability cache must really probe this older binary and detect no API.
capture macro drop rdistnn_distance_plugin_api
global rdistnn_distance_engine auto
mata: nc_got=rdistnn_distance_matrix(st_matrix("NC_X"),st_matrix("NC_T"),"l2");st_numscalar("nc_equal",max(abs(nc_ref-nc_got)))
test_assert "old plugin capability probe marks unsupported version" "$rdistnn_distance_plugin_api"=="0"
test_assert "old plugin auto falls back with reference distances" "$rdistnn_last_distance_engine"=="mata" & scalar(nc_equal)==0
global rdistnn_use_plugin 1
capture mata: nc_got=rdistnn_distance_matrix(st_matrix("NC_X"),st_matrix("NC_T"),"l2")
test_assert "old plugin auto with shared forcing fails explicitly" _rc!=0
capture macro drop rdistnn_use_plugin
global rdistnn_distance_engine plugin
capture mata: nc_got=rdistnn_distance_matrix(st_matrix("NC_X"),st_matrix("NC_T"),"l2")
test_assert "old plugin forced route fails explicitly" _rc!=0
// Pretend the successful probe was cached before binary replacement. The real
// old plugin rejects distance_batch, so test actual operation failure as well.
global rdistnn_distance_plugin_api 1
global rdistnn_distance_engine auto
mata: nc_got=rdistnn_distance_matrix(st_matrix("NC_X"),st_matrix("NC_T"),"l2");st_numscalar("nc_equal",max(abs(nc_ref-nc_got)))
test_assert "native operation failure auto falls back" "$rdistnn_last_distance_engine"=="mata" & scalar(nc_equal)==0
global rdistnn_use_plugin 1
global rdistnn_distance_plugin_api 1
capture mata: nc_got=rdistnn_distance_matrix(st_matrix("NC_X"),st_matrix("NC_T"),"l2")
test_assert "native operation failure auto with shared forcing fails" _rc!=0
capture macro drop rdistnn_use_plugin
global rdistnn_distance_plugin_api 1
global rdistnn_distance_engine plugin
capture mata: nc_got=rdistnn_distance_matrix(st_matrix("NC_X"),st_matrix("NC_T"),"l2")
test_assert "native operation failure forced route fails explicitly" _rc!=0
capture macro drop rdistnn_distance_engine
capture macro drop rdistnn_last_distance_engine
capture macro drop rdistnn_jackknife_engine
capture macro drop rdistnn_use_plugin
capture macro drop rdistnn_distance_plugin_api
