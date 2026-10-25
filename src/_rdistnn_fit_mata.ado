*! version 0.1.5  2026-10-25
version 17
// Fresh-command IO adapters. Capture validated ado locals at entry, keep all
// numeric work local, and publish only the caller's temporary Stata outputs.
// Base kernels are loaded first; estimator extensions retain their own owners.
capture mata: mata drop rdistnn_np_fit_ado()
capture mata: mata drop rdistnn_cate_fit_ado()
capture mata: mata drop rdistnn_weights_fit_ado()
mata:
void rdistnn_np_fit_ado()
{
    real scalar ado_default_index,ado_k,ado_n,ado_ncandidates,ado_ntargets,ado_s,ado_scaleeps,_np_g,_np_h
    real scalar _np_prior,_np_sg
    string scalar ado_V,ado_b,ado_candidate_mat,ado_distance,ado_grid_diag,ado_grid_est,ado_grid_se
    string scalar ado_method,ado_model_distances,ado_model_geometry,ado_model_groups,ado_model_ranks
    string scalar ado_model_signals,ado_neff_mat,ado_outcome,ado_scalerule,ado_scalevec,ado_se
    string scalar ado_target_diag,ado_targets_mat,ado_touse,ado_xvars
    real matrix _np_DX,_np_V,_np_Vg,_np_W,_np_Wg,_np_X,_np_X0,_np_candidates,_np_geometry,_np_grid_diag,_np_grid_est,_np_grid_se,_np_groups,_np_mu,_np_neff,_np_ranks,_np_target_diag,_np_y
    ado_default_index=strtoreal(st_local("default_index"))
    ado_k=strtoreal(st_local("k"))
    ado_n=strtoreal(st_local("n"))
    ado_ncandidates=strtoreal(st_local("ncandidates"))
    ado_ntargets=strtoreal(st_local("ntargets"))
    ado_s=strtoreal(st_local("s"))
    ado_scaleeps=strtoreal(st_local("scaleeps"))
    ado_V=st_local("V")
    ado_b=st_local("b")
    ado_candidate_mat=st_local("candidate_mat")
    ado_distance=st_local("distance")
    ado_grid_diag=st_local("grid_diag")
    ado_grid_est=st_local("grid_est")
    ado_grid_se=st_local("grid_se")
    ado_method=st_local("method")
    ado_model_distances=st_local("model_distances")
    ado_model_geometry=st_local("model_geometry")
    ado_model_groups=st_local("model_groups")
    ado_model_ranks=st_local("model_ranks")
    ado_model_signals=st_local("model_signals")
    ado_neff_mat=st_local("neff_mat")
    ado_outcome=st_local("outcome")
    ado_scalerule=st_local("scalerule")
    ado_scalevec=st_local("scalevec")
    ado_se=st_local("se")
    ado_target_diag=st_local("target_diag")
    ado_targets_mat=st_local("targets_mat")
    ado_touse=st_local("touse")
    ado_xvars=st_local("xvars")

    _np_X = st_data(., ado_xvars, ado_touse)
    _np_X0 = st_matrix(ado_targets_mat)
    _np_geometry = _np_DX = _np_ranks = J(0,0,.)
    rdistnn_localization_geometry(_np_X,_np_X0,ado_distance,ado_scalerule,
        strtoreal(tokens(ado_scalevec)),ado_scaleeps,_np_geometry,_np_DX,_np_ranks)
    _np_groups=rdistnn_distance_groups(_np_ranks,_np_DX,1)
    _np_y = st_data(., ado_outcome, ado_touse)
    _np_candidates = st_matrix(ado_candidate_mat)
    _np_grid_est = J(ado_ncandidates, ado_ntargets, .)
    _np_grid_se = J(ado_ncandidates, ado_ntargets, .)
    _np_grid_diag = J(ado_ncandidates, 5, .)

    for (_np_g = 1; _np_g <= ado_ncandidates; _np_g++) {
        _np_grid_diag[_np_g, 3] = min((5, ado_n))
        if (_np_candidates[_np_g, 4] == 0) continue
        _np_prior = 0
        if (_np_g > 1) {
            for (_np_h = 1; _np_h < _np_g; _np_h++) {
                if (_np_candidates[_np_h, 4] == 1 &
                    _np_candidates[_np_h, 2] == _np_candidates[_np_g, 2]) {
                    _np_prior = _np_h
                    break
                }
            }
        }
        if (_np_prior > 0) {
            _np_grid_est[_np_g,.] = _np_grid_est[_np_prior,.]
            _np_grid_se[_np_g,.] = _np_grid_se[_np_prior,.]
            _np_grid_diag[_np_g,.] = _np_grid_diag[_np_prior,.]
        }
        else {
            _np_sg = _np_candidates[_np_g, 8]
            _np_Wg = rdistnn_weights_from_rmat(_np_ranks, ado_k, _np_sg, ado_method,_np_DX,_np_groups)
            _np_grid_est[_np_g,.] = rdistnn_weighted_signals(_np_y, _np_Wg)
            _np_grid_diag[_np_g,.]=rdistnn_grid_weight_diagnostics(_np_Wg)
            if (ado_se == "") {
                _np_Vg = rdistnn_jk_cov_multi_from_ranks(
                    _np_ranks, _np_y, ado_k, _np_sg, ado_method,_np_DX,_np_groups)
                _np_grid_se[_np_g,.] = sqrt(diagonal(_np_Vg))'
            }
            // Save the default's first representative, including when the
            // default row duplicates an earlier candidate after other orders.
            if (_np_candidates[_np_g, 2] == _np_candidates[ado_default_index, 2]) {
                _np_W = _np_Wg
                if (ado_se == "") _np_V = _np_Vg
            }
        }
    }

    _np_mu = _np_grid_est[ado_default_index,.]
    _np_neff = 1 :/ colsum(_np_W :^ 2)
    _np_target_diag = rdistnn_diagnostics("np_dnn",_np_X,_np_X0,
        _np_y,_np_W,_np_DX,ado_k,ado_s,_np_mu,J(0,1,.),J(0,1,.))
    st_matrix(ado_b, _np_mu)
    st_matrix(ado_neff_mat, _np_neff)
    st_matrix(ado_grid_est, _np_grid_est)
    st_matrix(ado_grid_se, _np_grid_se)
    st_matrix(ado_grid_diag, _np_grid_diag)
    st_matrix(ado_target_diag, _np_target_diag)
    st_numscalar("r(neff_min)", min(_np_neff))
    st_numscalar("r(neff_max)", max(_np_neff))
    st_matrix(ado_model_groups,_np_groups); st_matrix(ado_model_signals,_np_y); st_matrix(ado_model_ranks,_np_ranks); st_matrix(ado_model_distances,_np_DX)
    if(any(ado_distance:==("mahalanobis","l1_scaled","l2_scaled"))) st_matrix(ado_model_geometry,_np_geometry)
    if(ado_se=="") st_matrix(ado_V,_np_V)

}

void rdistnn_cate_fit_ado()
{
    real scalar ado_default_index,ado_k,ado_n,ado_ncandidates,ado_ntargets,ado_s,ado_scaleeps,_ct_g,_ct_h
    real scalar _ct_prior,_ct_sg
    string scalar ado_V,ado_b,ado_candidate_mat,ado_distance,ado_eclip,ado_grid_diag,ado_grid_est
    string scalar ado_grid_se,ado_method,ado_model_distances,ado_model_geometry,ado_model_groups
    string scalar ado_model_ranks,ado_model_signals,ado_neff_mat,ado_rho_use,ado_scalerule,ado_scalevec
    string scalar ado_score_method,ado_se,ado_target_diag,ado_targets_mat,ado_touse,ado_treat_var
    string scalar ado_xvars
    real matrix _ct_DX,_ct_E,_ct_T,_ct_V,_ct_Vg,_ct_W,_ct_Wg,_ct_X,_ct_X0,_ct_candidates,_ct_geometry,_ct_grid_diag,_ct_grid_est,_ct_grid_se,_ct_groups,_ct_neff,_ct_ranks,_ct_rho,_ct_target_diag,_ct_tau
    ado_default_index=strtoreal(st_local("default_index"))
    ado_k=strtoreal(st_local("k"))
    ado_n=strtoreal(st_local("n"))
    ado_ncandidates=strtoreal(st_local("ncandidates"))
    ado_ntargets=strtoreal(st_local("ntargets"))
    ado_s=strtoreal(st_local("s"))
    ado_scaleeps=strtoreal(st_local("scaleeps"))
    ado_V=st_local("V")
    ado_b=st_local("b")
    ado_candidate_mat=st_local("candidate_mat")
    ado_distance=st_local("distance")
    ado_eclip=st_local("eclip")
    ado_grid_diag=st_local("grid_diag")
    ado_grid_est=st_local("grid_est")
    ado_grid_se=st_local("grid_se")
    ado_method=st_local("method")
    ado_model_distances=st_local("model_distances")
    ado_model_geometry=st_local("model_geometry")
    ado_model_groups=st_local("model_groups")
    ado_model_ranks=st_local("model_ranks")
    ado_model_signals=st_local("model_signals")
    ado_neff_mat=st_local("neff_mat")
    ado_rho_use=st_local("rho_use")
    ado_scalerule=st_local("scalerule")
    ado_scalevec=st_local("scalevec")
    ado_score_method=st_local("score_method")
    ado_se=st_local("se")
    ado_target_diag=st_local("target_diag")
    ado_targets_mat=st_local("targets_mat")
    ado_touse=st_local("touse")
    ado_treat_var=st_local("treat_var")
    ado_xvars=st_local("xvars")

    _ct_X     = st_data(., ado_xvars, ado_touse)
    _ct_X0    = st_matrix(ado_targets_mat)
    _ct_geometry = _ct_DX = _ct_ranks = J(0,0,.)
    rdistnn_localization_geometry(_ct_X,_ct_X0,ado_distance,ado_scalerule,
        strtoreal(tokens(ado_scalevec)),ado_scaleeps,_ct_geometry,_ct_DX,_ct_ranks)
    _ct_groups=rdistnn_distance_groups(_ct_ranks,_ct_DX,1)
    _ct_rho   = st_data(., ado_rho_use, ado_touse)
    _ct_candidates = st_matrix(ado_candidate_mat)
    _ct_grid_est = J(ado_ncandidates, ado_ntargets, .)
    _ct_grid_se = J(ado_ncandidates, ado_ntargets, .)
    _ct_grid_diag = J(ado_ncandidates, 5, .)

    for (_ct_g = 1; _ct_g <= ado_ncandidates; _ct_g++) {
        _ct_grid_diag[_ct_g, 3] = min((5, ado_n))
        if (_ct_candidates[_ct_g, 4] == 0) continue
        _ct_prior = 0
        if (_ct_g > 1) {
            for (_ct_h = 1; _ct_h < _ct_g; _ct_h++) {
                if (_ct_candidates[_ct_h, 4] == 1 &
                    _ct_candidates[_ct_h, 2] == _ct_candidates[_ct_g, 2]) {
                    _ct_prior = _ct_h
                    break
                }
            }
        }
        if (_ct_prior > 0) {
            _ct_grid_est[_ct_g,.] = _ct_grid_est[_ct_prior,.]
            _ct_grid_se[_ct_g,.] = _ct_grid_se[_ct_prior,.]
            _ct_grid_diag[_ct_g,.] = _ct_grid_diag[_ct_prior,.]
        }
        else {
            _ct_sg = _ct_candidates[_ct_g, 8]
            _ct_Wg = rdistnn_weights_from_rmat(_ct_ranks, ado_k, _ct_sg, ado_method,_ct_DX,_ct_groups)
            _ct_grid_est[_ct_g,.] = rdistnn_weighted_signals(_ct_rho, _ct_Wg)
            _ct_grid_diag[_ct_g,.]=rdistnn_grid_weight_diagnostics(_ct_Wg)
            if (ado_se == "") {
                _ct_Vg = rdistnn_jk_cov_multi_from_ranks(
                    _ct_ranks, _ct_rho, ado_k, _ct_sg, ado_method,_ct_DX,_ct_groups)
                _ct_grid_se[_ct_g,.] = sqrt(diagonal(_ct_Vg))'
            }
            // Save the default's first representative, including when the
            // default row duplicates an earlier candidate after other orders.
            if (_ct_candidates[_ct_g, 2] == _ct_candidates[ado_default_index, 2]) {
                _ct_W = _ct_Wg
                if (ado_se == "") _ct_V = _ct_Vg
            }
        }
    }

    _ct_tau = _ct_grid_est[ado_default_index,.]
    _ct_neff = 1 :/ colsum(_ct_W :^ 2)
    _ct_T = _ct_E = J(0,1,.)
    if (ado_score_method != "pseudo_outcome") {
        _ct_T = st_data(., ado_treat_var, ado_touse)
        _ct_E = st_data(., ado_eclip, ado_touse)
    }
    _ct_target_diag = rdistnn_diagnostics("cate_dnn",_ct_X,_ct_X0,
        _ct_rho,_ct_W,_ct_DX,ado_k,ado_s,_ct_tau,_ct_T,_ct_E)
    st_matrix(ado_b, _ct_tau)
    st_matrix(ado_neff_mat, _ct_neff)
    st_matrix(ado_grid_est, _ct_grid_est)
    st_matrix(ado_grid_se, _ct_grid_se)
    st_matrix(ado_grid_diag, _ct_grid_diag)
    st_matrix(ado_target_diag, _ct_target_diag)
    st_numscalar("r(neff_min)", min(_ct_neff))
    st_numscalar("r(neff_max)", max(_ct_neff))
    st_matrix(ado_model_groups,_ct_groups); st_matrix(ado_model_signals,_ct_rho); st_matrix(ado_model_ranks,_ct_ranks); st_matrix(ado_model_distances,_ct_DX)
    if(any(ado_distance:==("mahalanobis","l1_scaled","l2_scaled"))) st_matrix(ado_model_geometry,_ct_geometry)
    if(ado_se=="") st_matrix(ado_V,_ct_V)

}

void rdistnn_weights_fit_ado()
{
    real scalar ado_k,ado_s,ado_scaleeps
    string scalar ado_distance,ado_distances,ado_geometry,ado_groups,ado_method,ado_neff,ado_ranks
    string scalar ado_row_ids,ado_scalerule,ado_scalevec,ado_targets,ado_touse,ado_varlist,ado_weights
    real matrix _wd_D,_wd_G,_wd_R,_wd_W,_wd_X,_wd_X0,_wd_groups
    ado_k=strtoreal(st_local("k"))
    ado_s=strtoreal(st_local("s"))
    ado_scaleeps=strtoreal(st_local("scaleeps"))
    ado_distance=st_local("distance")
    ado_distances=st_local("distances")
    ado_geometry=st_local("geometry")
    ado_groups=st_local("groups")
    ado_method=st_local("method")
    ado_neff=st_local("neff")
    ado_ranks=st_local("ranks")
    ado_row_ids=st_local("row_ids")
    ado_scalerule=st_local("scalerule")
    ado_scalevec=st_local("scalevec")
    ado_targets=st_local("targets")
    ado_touse=st_local("touse")
    ado_varlist=st_local("varlist")
    ado_weights=st_local("weights")

    _wd_X=st_data(.,tokens(ado_varlist),ado_touse)
    _wd_X0=st_matrix(ado_targets)
    _wd_G=rdistnn_fit_geometry(_wd_X,ado_distance,ado_scalerule,strtoreal(tokens(ado_scalevec)),ado_scaleeps)
    _wd_D=rdistnn_distance_matrix(_wd_X,_wd_X0,ado_distance,_wd_G)
    _wd_R=rdistnn_ranks_from_dmat(_wd_D)
    _wd_groups=rdistnn_distance_groups(_wd_R,_wd_D,1)
    rdistnn_distance_warning(_wd_groups)
    _wd_W=rdistnn_weights_from_rmat(_wd_R,ado_k,ado_s,ado_method,_wd_D,_wd_groups)
    st_matrix(ado_weights,_wd_W); st_matrix(ado_ranks,_wd_R)
    st_matrix(ado_distances,_wd_D)
    st_matrix(ado_groups,_wd_groups)
    st_matrix(ado_row_ids,select((1::st_nobs()),st_data(.,ado_touse)))
    st_matrix(ado_neff,1:/colsum(_wd_W:^2))
    if(any(ado_distance:==("mahalanobis","l1_scaled","l2_scaled"))) st_matrix(ado_geometry,_wd_G)

}
end
