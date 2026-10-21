version 17
// Frozen geometry, ordered distances and ranks. Estimator scores remain separate.
capture mata: mata drop rdistnn_localization_geometry()
capture mata: mata drop rdistnn_localization_ranks()
mata:
// One covariate-aware rank boundary for fresh, restored and plotting routes.
real matrix rdistnn_localization_ranks(real matrix DX, real matrix X)
{
    return(rdistnn_ranks_from_dmat(DX,X))
}

void rdistnn_localization_geometry(real matrix X, real matrix X0,
    string scalar distance, string scalar scalerule, real rowvector scalevec,
    real scalar scaleeps, real matrix geometry, real matrix DX, real matrix ranks)
{
    geometry=rdistnn_fit_geometry(X,distance,scalerule,scalevec,scaleeps)
    DX=rdistnn_distance_matrix(X,X0,distance,geometry)
    ranks=rdistnn_localization_ranks(DX,X)
}
end
