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
    capture program define rdistnn_plugin, plugin using(rdistnn.plugin)
    local plugin_rc=_rc
    // At this declaration boundary, 601 means the file is absent and 9999
    // means the dynamic loader could not load it (for example, a host mismatch).
    // Other failures, including Break and resource errors, are not capability.
    if inlist(`plugin_rc',0,110) {
        global rdistnn_plugin_ok 1
    }
    else if inlist(`plugin_rc',601,9999) {
        global rdistnn_plugin_ok 0
    }
    else {
        display as error "Native plugin declaration failed (return code `plugin_rc'). Resolve the execution or resource error and rerun the command."
        exit `plugin_rc'
    }
}

foreach module in messages native_transport accumulation covariance blocks geometry distance_plugin distance_mata rank_weights rank_loo localization_mata diagnostic_mata plqr_mata plot_mata fit_mata {
    findfile "_rdistnn_`module'.ado"
    quietly run "`r(fn)'"
}
