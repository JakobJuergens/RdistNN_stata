version 17
// Shared gateway for the Mata core. Load dependencies in order;
// each module owns its functions and reload guards. Estimator commands continue
// to source this single entry point in fresh and existing sessions.
//
// To build an mlib: run this file, cd src, then use
// mata: mata mlib create lrdistnn, replace
// mata: mata mlib add lrdistnn rdistnn_*()

// Declare the optional plugin without changing caller-selected engines.
if ("$rdistnn_plugin_ok" == "") {
    // _rc == 0   : declared for the first time
    // _rc == 110 : already declared in this session
    // anything else: rdistnn.plugin file not found on adopath
    capture program define rdistnn_plugin, plugin using(rdistnn.plugin)
    if _rc == 0 | _rc == 110 {
        global rdistnn_plugin_ok 1
    }
    else {
        global rdistnn_plugin_ok 0
    }
}

foreach module in native_transport accumulation covariance blocks geometry distance_plugin distance_mata rank_weights rank_loo localization_mata diagnostic_mata plqr_mata plot_mata {
    findfile "_rdistnn_`module'.ado"
    quietly run "`r(fn)'"
}
