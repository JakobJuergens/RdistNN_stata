*! version 0.1.5  2026-10-25
version 17
// Paired CLATE evaluator. The caller loads _rdistnn_mata first; reloading
// this extension drops only its own functions and preserves the base library.
capture mata: mata drop rdistnn_clate_fit()

mata:

// Last five arguments are outputs. componentV is interleaved Y1,D1,Y2,D2.
// Predictions, ranks, metric and s remain fixed across paired deletions.
void rdistnn_clate_fit(real matrix ranks, real matrix signals,
    real scalar k, real scalar s, string scalar method,
    real scalar do_se, real scalar tol, real matrix components,
    real rowvector theta, real matrix V, real matrix componentV,
    real matrix status, | real matrix DX, real matrix prepared_groups)
{
    real scalar n, m, j, h, d, den_scale, den_rescale, scale, variance, corr, omega_quarter
    real colvector sd, good
    real matrix W, loo, centered, g, normalized, block

    if(args()<14) prepared_groups=J(0,0,.)
    n=rows(ranks); m=cols(ranks)
    if (n<1 | m<1 | rows(signals)!=n | cols(signals)!=2 |
        hasmissing(ranks) | hasmissing(signals) | tol<0 | tol>=. |
        k<1 | k!=floor(k) | s<k | s!=floor(s) | s>n |
        (do_se & s>=n)) {
        rdistnn_fail(198,"rdistnn_clate_fit: invalid ranks, signals, tuning, or tolerance.")
    }
    if(args()>=13 & rows(DX)>0) W=rdistnn_weights_from_rmat(ranks,k,s,method,DX,prepared_groups)
    else W=rdistnn_weights_from_rmat(ranks,k,s,method)
    components=rdistnn_weighted_signals(signals,W)'
    theta=J(1,m,.)
    V=J(m,m,.)
    componentV=J(2*m,2*m,.)
    status=J(m,5,0)
    status[,4]=J(m,1,1)
    // Mata's negative nonmissing range is twice its positive range. Detect
    // unsafe magnitudes before abs(), whose missing result max() would ignore.
    // Halving both sides is exact binary scaling and retains inactive signals.
    den_rescale=(min(signals[,2]) < -maxdouble() |
        min(components[,2]) < -maxdouble())
    if(den_rescale) den_scale=max(abs(signals[,2]/2))
    else den_scale=max(abs(signals[,2]))
    for(j=1;j<=m;j++) {
        d=components[j,2]
        status[j,2]=(d<0)
        if(d!=0 & d<. & den_scale>0 & tol>0) {
            if(den_rescale) status[j,3]=(abs(d/2)/den_scale<=tol)
            else status[j,3]=(abs(d)/den_scale<=tol)
        }
        if(hasmissing(components[j,])) status[j,1]=2
        else if(d==0) status[j,1]=1
        else {
            theta[j]=components[j,1]/d
            if(theta[j]>=.) status[j,1]=2
        }
    }
    if(!do_se) return

    if(args()>=13 & rows(DX)>0) {
        loo=rdistnn_loo_from_ranks(ranks,signals,k,s,method,1,
            rdistnn_weights_from_ranks((1::(n-1)),k,s,method),DX,prepared_groups)
    }
    else loo=rdistnn_loo_from_ranks(ranks,signals,k,s,method,1)
    centered=loo:-mean(loo)
    componentV=rdistnn_covariance(loo,(n-1)/n,1,1)
    normalized=J(n,m,.)
    sd=J(m,1,.)
    good=J(m,1,0)
    status[,4]=J(m,1,0)
    for(j=1;j<=m;j++) {
        block=centered[,(2*j-1)..(2*j)]
        if(hasmissing(block) |
            hasmissing(componentV[(2*j-1)..(2*j),(2*j-1)..(2*j)]) |
            (any(block[,1]:!=0) & componentV[2*j-1,2*j-1]==0) |
            (any(block[,2]:!=0) & componentV[2*j,2*j]==0)) {
            componentV[(2*j-1)..(2*j),]=J(2,2*m,.)
            componentV[,(2*j-1)..(2*j)]=J(2*m,2,.)
            status[j,4]=3
            continue
        }
        if(theta[j]>=.) continue
        // Linearize centered component deletions at the full-sample ratio.
        // Divide before combining to avoid squaring a tiny denominator.
        g=block[,1]/components[j,2]-theta[j]*(block[,2]/components[j,2])
        scale=max(abs(g))
        if(hasmissing(g) | scale>=.) {
            status[j,4]=3
            continue
        }
        if(scale==0) {
            normalized[,j]=J(n,1,0)
            sd[j]=0
            status[j,4]=2
        }
        else {
            normalized[,j]=g/scale
            sd[j]=sqrt(rdistnn_covariance(g,(n-1)/n,0,1)[1,1])
        }
        variance=sd[j]^2
        if(variance>=. | (scale>0 & variance==0)) {
            status[j,4]=3
            continue
        }
        // R also requires finite neff*variance. Stata's nonmissing range is
        // smaller than R's IEEE range, so compare one quarter of that value
        // with DBL_MAX/4 instead of forming an unsupported omega directly.
        omega_quarter=variance*((1/quadsum(W[,j]:^2))/4)
        if(omega_quarter>4.4942328371557893e307) {
            status[j,4]=3
            continue
        }
        V[j,j]=variance
        good[j]=1
        status[j,5]=(variance>0)
    }
    for(j=1;j<=m;j++) if(good[j]) {
        for(h=j+1;h<=m;h++) if(good[h]) {
            if(sd[j]==0 | sd[h]==0) V[j,h]=0
            else {
                corr=cross(normalized[,j],normalized[,h]) /
                    sqrt(quadsum(normalized[,j]:^2)*quadsum(normalized[,h]:^2))
                V[j,h]=sd[j]*(corr*sd[h])
            }
            V[h,j]=V[j,h]
        }
    }
}


end

// Fresh-command IO with function-local scratch.
capture mata: mata drop rdistnn_clate_fit_ado()
mata:
void rdistnn_clate_fit_ado()
{
    real scalar ado_active,ado_default,ado_denominatortol,ado_k,ado_m,ado_n,ado_ng,ado_resolved_s
    real scalar ado_scaleeps,_cl_g,_cl_h,_cl_prior
    string scalar ado_acol,ado_componentV,ado_components,ado_diag,ado_eclip,ado_grid_diag,ado_gridd
    string scalar ado_gridse,ado_gridstatus,ado_gridtheta,ado_gridy,ado_meta,ado_model_distances
    string scalar ado_model_geometry,ado_model_groups,ado_model_ranks,ado_model_signals,ado_neff,ado_phid
    string scalar ado_phiy,ado_resolved_distance,ado_resolved_method,ado_route,ado_scalerule,ado_scalevec
    string scalar ado_semat,ado_status,ado_targets,ado_theta,ado_touse,ado_vcov,ado_xvars
    real rowvector _cl_diag_representative
    real matrix _cl_A,_cl_CV,_cl_DX,_cl_V,_cl_W,_cl_X,_cl_X0,_cl_b,_cl_c,_cl_e,_cl_geometry,_cl_grid,_cl_grid_W,_cl_grid_diag,_cl_gridd,_cl_gridse,_cl_gridstatus,_cl_gridy,_cl_groups,_cl_meta,_cl_ranks,_cl_signals,_cl_status
    ado_active=strtoreal(st_local("active"))
    ado_default=strtoreal(st_local("default"))
    ado_denominatortol=strtoreal(st_local("denominatortol"))
    ado_k=strtoreal(st_local("k"))
    ado_m=strtoreal(st_local("m"))
    ado_n=strtoreal(st_local("n"))
    ado_ng=strtoreal(st_local("ng"))
    ado_resolved_s=strtoreal(st_local("resolved_s"))
    ado_scaleeps=strtoreal(st_local("scaleeps"))
    ado_acol=st_local("acol")
    ado_componentV=st_local("componentV")
    ado_components=st_local("components")
    ado_diag=st_local("diag")
    ado_eclip=st_local("eclip")
    ado_grid_diag=st_local("grid_diag")
    ado_gridd=st_local("gridd")
    ado_gridse=st_local("gridse")
    ado_gridstatus=st_local("gridstatus")
    ado_gridtheta=st_local("gridtheta")
    ado_gridy=st_local("gridy")
    ado_meta=st_local("meta")
    ado_model_distances=st_local("model_distances")
    ado_model_geometry=st_local("model_geometry")
    ado_model_groups=st_local("model_groups")
    ado_model_ranks=st_local("model_ranks")
    ado_model_signals=st_local("model_signals")
    ado_neff=st_local("neff")
    ado_phid=st_local("phid")
    ado_phiy=st_local("phiy")
    ado_resolved_distance=st_local("resolved_distance")
    ado_resolved_method=st_local("resolved_method")
    ado_route=st_local("route")
    ado_scalerule=st_local("scalerule")
    ado_scalevec=st_local("scalevec")
    ado_semat=st_local("semat")
    ado_status=st_local("status")
    ado_targets=st_local("targets")
    ado_theta=st_local("theta")
    ado_touse=st_local("touse")
    ado_vcov=st_local("vcov")
    ado_xvars=st_local("xvars")

    _cl_X = st_data(., ado_xvars, ado_touse)
    _cl_X0 = st_matrix(ado_targets)
    _cl_geometry = _cl_DX = _cl_ranks = J(0,0,.)
    rdistnn_localization_geometry(_cl_X,_cl_X0,ado_resolved_distance,ado_scalerule,
        strtoreal(tokens(ado_scalevec)),ado_scaleeps,_cl_geometry,_cl_DX,_cl_ranks)
    _cl_groups=rdistnn_distance_groups(_cl_ranks,_cl_DX,1)
    _cl_signals = st_data(., (ado_phiy+" "+ado_phid), ado_touse)
    _cl_meta = st_matrix(ado_meta)
    _cl_grid = _cl_gridse = _cl_gridy = _cl_gridd = J(ado_ng,ado_m,.)
    _cl_gridstatus = J(ado_ng*ado_m,7,.)
    _cl_c = _cl_V = _cl_CV = _cl_status = J(0,0,.)
    _cl_b = J(1,0,.)
    for (_cl_g=1; _cl_g<=ado_ng; _cl_g++) {
        if (_cl_meta[_cl_g,4]==0) continue
        _cl_prior=0
        for (_cl_h=1; _cl_h<_cl_g; _cl_h++) {
            if (_cl_meta[_cl_h,4] & _cl_meta[_cl_h,2]==_cl_meta[_cl_g,2]) {
                _cl_prior=_cl_h
                break
            }
        }
        if (_cl_prior) {
            _cl_grid[_cl_g,.]=_cl_grid[_cl_prior,.]
            _cl_gridse[_cl_g,.]=_cl_gridse[_cl_prior,.]
            _cl_gridy[_cl_g,.]=_cl_gridy[_cl_prior,.]
            _cl_gridd[_cl_g,.]=_cl_gridd[_cl_prior,.]
            _cl_gridstatus[((_cl_g-1)*ado_m+1)..(_cl_g*ado_m),.] =
                (J(ado_m,1,_cl_g),_cl_gridstatus[((_cl_prior-1)*ado_m+1)..(_cl_prior*ado_m),2..7])
            continue
        }
        rdistnn_clate_fit(_cl_ranks,_cl_signals,ado_k,_cl_meta[_cl_g,8],
            ado_resolved_method,ado_active,ado_denominatortol,
            _cl_c,_cl_b,_cl_V,_cl_CV,_cl_status,_cl_DX,_cl_groups)
        _cl_grid[_cl_g,.] = _cl_b
        _cl_gridse[_cl_g,.] = sqrt(diagonal(_cl_V))'
        _cl_gridy[_cl_g,.] = _cl_c[,1]'
        _cl_gridd[_cl_g,.] = _cl_c[,2]'
        _cl_gridstatus[((_cl_g-1)*ado_m+1)..(_cl_g*ado_m),.] =
            (J(ado_m,1,_cl_g),(1::ado_m),_cl_status)
        if (_cl_meta[_cl_g,2]==_cl_meta[ado_default,2]) {
            st_matrix(ado_components,_cl_c)
            st_matrix(ado_theta,_cl_b)
            st_matrix(ado_vcov,_cl_V)
            st_matrix(ado_componentV,_cl_CV)
            st_matrix(ado_status,_cl_status)
        }
    }
    _cl_b = st_matrix(ado_theta)
    _cl_W = rdistnn_weights_from_rmat(_cl_ranks,ado_k,ado_resolved_s,ado_resolved_method,_cl_DX,_cl_groups)
    _cl_A = _cl_e = J(0,1,.)
    if (ado_route=="augmented_iv") {
        _cl_A = st_data(.,ado_acol,ado_touse)
        _cl_e = st_data(.,ado_eclip,ado_touse)
    }
    st_matrix(ado_diag,rdistnn_diagnostics("clate_dnn",_cl_X,_cl_X0,
        _cl_signals,_cl_W,_cl_DX,ado_k,ado_resolved_s,_cl_b,_cl_A,_cl_e))
    st_matrix(ado_neff,1:/colsum(_cl_W:^2))
    _cl_grid_diag=J(ado_ng,5,.)
    _cl_diag_representative=J(1,ado_ng,0)
    for(_cl_g=1;_cl_g<=ado_ng;_cl_g++) {
        _cl_grid_diag[_cl_g,3]=min((5,ado_n))
        if(_cl_meta[_cl_g,4]==0) continue
        _cl_h=_cl_meta[_cl_g,2]
        _cl_prior=_cl_diag_representative[_cl_h]
        if(_cl_prior) _cl_grid_diag[_cl_g,.]=_cl_grid_diag[_cl_prior,.]
        else {
            // Aliases share diagnostics and the default already has its W.
            // Retain only one scratch matrix for other distinct orders.
            if(_cl_h==_cl_meta[ado_default,2]) {
                _cl_grid_diag[_cl_g,.]=rdistnn_grid_weight_diagnostics(_cl_W)
            }
            else {
                _cl_grid_W=rdistnn_weights_from_rmat(_cl_ranks,ado_k,_cl_meta[_cl_g,8],ado_resolved_method,_cl_DX,_cl_groups)
                _cl_grid_diag[_cl_g,.]=rdistnn_grid_weight_diagnostics(_cl_grid_W)
            }
            _cl_diag_representative[_cl_h]=_cl_g
        }
    }
    st_matrix(ado_grid_diag,_cl_grid_diag)
    st_matrix(ado_semat,sqrt(diagonal(st_matrix(ado_vcov)))')
    st_matrix(ado_gridtheta,_cl_grid)
    st_matrix(ado_gridse,_cl_gridse)
    st_matrix(ado_gridy,_cl_gridy)
    st_matrix(ado_gridd,_cl_gridd)
    st_matrix(ado_gridstatus,_cl_gridstatus)
    st_matrix(ado_model_groups,_cl_groups); st_matrix(ado_model_signals,_cl_signals); st_matrix(ado_model_ranks,_cl_ranks); st_matrix(ado_model_distances,_cl_DX)
    if(any(ado_resolved_distance:==("mahalanobis","l1_scaled","l2_scaled"))) st_matrix(ado_model_geometry,_cl_geometry)

}
end

capture mata: mata drop rdistnn_clate_post_ado()
mata:
void rdistnn_clate_post_ado()
{
    real rowvector keep
    real matrix covariance
    keep=strtoreal(tokens(st_local("keep")))
    st_matrix(st_local("postb"),st_matrix(st_local("theta"))[,keep])
    covariance=st_matrix(st_local("vcov"))[keep,keep]
    st_matrix(st_local("postV"),covariance)
    st_local("postcov",strofreal(!hasmissing(covariance)))
}
end
