version 17
// Validate the current numeric snapshot before any cached consumer reads it.
// Cold report/diagnostic calls may have neither a dataset nor a Mata library.
program define _rdistnn_model_check
    version 17
    if e(model_version)!=4 | "`e(tie_rule)'"!="exact_distance_block_averaging_v1" {
        di as error "Saved fits require the current exact-distance averaging policy; refit with the current package."
        exit 498
    }
    tempname core_ready
    mata: st_numscalar("`core_ready'",findexternal("rdistnn_distance_groups()")!=NULL)
    if !scalar(`core_ready') {
        findfile "_rdistnn_mata.ado"
        quietly run "`r(fn)'"
    }
    findfile "_rdistnn_model_decode.ado"
    quietly run "`r(fn)'"
    mata: rdistnn_model_validate()
end
