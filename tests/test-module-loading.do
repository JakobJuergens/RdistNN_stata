// Cold loading and repeated loading must expose one complete shared core.
version 17
do "tests/test_framework.do"

test_section "shared modules: cold loading, reloads and independent references"
local old_distance "$rdistnn_distance_engine"
local old_jackknife "$rdistnn_jackknife_engine"
local old_plugin "$rdistnn_use_plugin"
global rdistnn_distance_engine mata
global rdistnn_jackknife_engine mata
global rdistnn_use_plugin 0

mata: mata clear
mata:
real scalar maintainability_foreign()
{
    return(123)
}
end

forvalues reload=1/2 {
    quietly do "src/_rdistnn_mata.ado"
    mata {
    ML_X = (0\1\3\6)
    ML_X0 = (0\4)
    ML_Y = (2\-1\4\7)
    ML_G = rdistnn_fit_geometry(ML_X,"l2","sd",J(1,0,.),1e-8)
    ML_D = rdistnn_distance_matrix(ML_X,ML_X0,"l2",ML_G)
    ML_R = rdistnn_ranks_from_dmat(ML_D)
    ML_W = rdistnn_weights_from_rmat(ML_R,1,2,"exact")
    ML_L = rdistnn_loo_from_ranks(ML_R,ML_Y,1,2,"exact")
    ML_V = rdistnn_jk_cov_multi_from_ranks(ML_R,ML_Y,1,2,"exact")
    st_numscalar("ml_distance",max(abs(ML_D-(0,4\1,3\3,1\6,2))))
    st_numscalar("ml_ranks",max(abs(ML_R-(1,4\2,3\3,1\4,2))))
    st_numscalar("ml_weights",max(abs(ML_W-(.5,0\1/3,1/6\1/6,.5\0,1/3))))
    st_numscalar("ml_estimate",max(abs(rdistnn_weighted_signals(ML_Y,ML_W)-(4/3,25/6))))
    st_numscalar("ml_deletions",max(abs(ML_L-(2/3,5\8/3,5\1,13/3\1,7/3))))
    st_numscalar("ml_covariance",max(abs(ML_V-(11/6,5/6\5/6,43/12))))
    st_numscalar("ml_quantile",rdistnn_weighted_quantile(ML_X,(0\.4\.6\0),.2))
    st_numscalar("ml_endpoint",rdistnn_weighted_quantile(ML_X,(0\.4\.6\0),0))
    st_numscalar("ml_foreign",maintainability_foreign())
    st_numscalar("ml_power",rdistnn_power2(-20))
    }
    test_assert "[modules `reload'] distances match original coordinates" scalar(ml_distance)==0
    test_assert "[modules `reload'] ranks preserve target columns" scalar(ml_ranks)==0
    test_assert "[modules `reload'] exact weights match the combinatorial law" scalar(ml_weights)<1e-14
    test_assert "[modules `reload'] weighted signals match hand calculation" scalar(ml_estimate)<1e-14
    test_assert "[modules `reload'] deletions match explicit subset means" scalar(ml_deletions)<1e-14
    test_assert "[modules `reload'] joint covariance matches explicit deletions" scalar(ml_covariance)<1e-14
    test_assert "[modules `reload'] positive quantile skips zero-weight support" scalar(ml_quantile)==1
    test_assert "[modules `reload'] zero probability keeps the first sorted observation" scalar(ml_endpoint)==0
    test_assert "[modules `reload'] unrelated caller function survives" scalar(ml_foreign)==123
    test_assert "[modules `reload'] shared binary scale remains available" scalar(ml_power)==2^(-20)
}

mata: mata drop maintainability_foreign()
mata: mata drop ML_X ML_X0 ML_Y ML_G ML_D ML_R ML_W ML_L ML_V
scalar drop ml_distance ml_ranks ml_weights ml_estimate ml_deletions ///
    ml_covariance ml_quantile ml_endpoint ml_foreign ml_power
global rdistnn_distance_engine "`old_distance'"
global rdistnn_jackknife_engine "`old_jackknife'"
global rdistnn_use_plugin "`old_plugin'"
