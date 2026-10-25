// Actual preparation counts and cached-target reuse through public commands.
// The fixture instruments production grouping, restores its modules/search path,
// and distinguishes retained targets from geometry that must be recomputed.
version 17
do "tests/test_framework.do"
test_section "group preparation is reused by public grids and cached consumers"
do "tests/helpers/reporting-instrument.do"
capture program drop _MR_group_reuse_fixture
program define _MR_group_reuse_fixture
    tempfile root
    local MR_SHADOW "`root'_sources"
    mkdir "`MR_SHADOW'"
    tempname base_fit base_cov base_diag base_grid_est base_grid_se
    capture noisily {
        mata: rdistnn_test_instrument("src/_rdistnn_blocks.ado",st_local("MR_SHADOW")+"/group_count.ado","n=rows(ranks); m=cols(ranks)","st_numscalar("+char(34)+"MR_GROUP_COUNT"+char(34)+",st_numscalar("+char(34)+"MR_GROUP_COUNT"+char(34)+")+1); n=rows(ranks); m=cols(ranks)")
        mata: rdistnn_test_instrument(st_local("MR_SHADOW")+"/group_count.ado",st_local("MR_SHADOW")+"/_rdistnn_blocks.ado","sorted=X[order", "st_numscalar("+char(34)+"MR_DUP_COUNT"+char(34)+",st_numscalar("+char(34)+"MR_DUP_COUNT"+char(34)+")+1); sorted=X[order")
        local weight_counter `"st_numscalar("MR_WEIGHT_COUNT",st_numscalar("MR_WEIGHT_COUNT")+1); if(s==4) st_numscalar("MR_DEFAULT_WEIGHT_COUNT",st_numscalar("MR_DEFAULT_WEIGHT_COUNT")+1); m = cols(ranks)"'
        mata: rdistnn_test_instrument("src/_rdistnn_rank_weights.ado",st_local("MR_SHADOW")+"/_rdistnn_rank_weights.ado","m = cols(ranks)",st_local("weight_counter"))
        adopath ++ "`MR_SHADOW'"

        clear
        set obs 24
        gen double x=mod(_n,6)-3
        gen double y=_n^2
        gen double den=1+_n/20
        matrix MR_TARGETS=(0\1\0)
        foreach command in np_dnn cate_dnn clate_dnn {
            local opts ""
            if "`command'"=="cate_dnn" local opts "pseudooutcome"
            if "`command'"=="clate_dnn" local opts "pseudooutcome denominator(den)"
            scalar MR_GROUP_COUNT=0
            scalar MR_DUP_COUNT=0
            scalar MR_WEIGHT_COUNT=0
            scalar MR_DEFAULT_WEIGHT_COUNT=0
            quietly `command' y, xcols(x) x0matrix(MR_TARGETS) tuningmode(grid) sgrid(3 4 5) sdefault(4) `opts'
            local distinct_weight_calls=scalar(MR_WEIGHT_COUNT)
            test_assert "`command': actual weight primitive counter has a positive control" `distinct_weight_calls'>0 & MR_DEFAULT_WEIGHT_COUNT>0
            test_assert "`command': default weights are reused by grid diagnostics" MR_DEFAULT_WEIGHT_COUNT<=2
            matrix `base_fit'=e(estimates)
            matrix `base_cov'=e(estimate_vcov)
            matrix `base_diag'=e(tuning_diagnostics)
            matrix `base_grid_est'=e(tuning_estimates)
            matrix `base_grid_se'=e(tuning_se)
            scalar MR_GROUP_COUNT=0
            scalar MR_DUP_COUNT=0
            scalar MR_WEIGHT_COUNT=0
            scalar MR_DEFAULT_WEIGHT_COUNT=0
            quietly `command' y, xcols(x) x0matrix(MR_TARGETS) tuningmode(grid) sgrid(3 4 3 5) sdefault(4) `opts'
            test_assert "`command': an aliased candidate adds no weight construction" MR_WEIGHT_COUNT==`distinct_weight_calls'
            mata: st_numscalar("MR_GRID_PUBLIC_SAME",all(st_matrix("e(estimates)"):==st_matrix(st_local("base_fit"))) & all(st_matrix("e(estimate_vcov)"):==st_matrix(st_local("base_cov"))) & all(st_matrix("e(tuning_diagnostics)"):==st_matrix(st_local("base_diag"))[(1\2\1\3),]) & all(st_matrix("e(tuning_estimates)"):==st_matrix(st_local("base_grid_est"))[(1\2\1\3),]) & all(st_matrix("e(tuning_se)"):==st_matrix(st_local("base_grid_se"))[(1\2\1\3),]) & all(st_matrix("e(tuning_candidates)")[,1]:==(1::4)) & all(st_matrix("e(tuning_candidates)")[,2]:==(1\2\1\3)))
            test_assert "`command': aliases preserve default outputs and every public candidate row" MR_GRID_PUBLIC_SAME==1
            test_assert "`command': one preparation across grid aliases kernels and posting" MR_GROUP_COUNT==1
            test_assert "`command': prepared fitting performs no duplicate-detection sorts" MR_DUP_COUNT==0
            scalar MR_GROUP_COUNT=0
            quietly rdistnn_predict, x0matrix(MR_TARGETS) sefit
            test_assert "`command': identical saved geometry reuses validated groups" MR_GROUP_COUNT==0
            scalar MR_GROUP_COUNT=0
            quietly rdistnn_predict, x0(0.5) sefit
            test_assert "`command': new geometry prepares its own groups once" MR_GROUP_COUNT==1
        }
        gen double mr_d=1
        gen double mr_zero=0
        foreach command in np_tdnn plqr_dnn {
            scalar MR_GROUP_COUNT=0
            scalar MR_DUP_COUNT=0
            if "`command'"=="np_tdnn" quietly np_tdnn y, xcols(x) x0matrix(MR_TARGETS) s1(2) s2(4)
            if "`command'"=="plqr_dnn" quietly plqr_dnn y, x(x) x0(0) d(mr_d) g(mr_zero) md(mr_zero) thetaeval(0) s(4)
            test_assert "`command': weights and covariance reuse one prepared geometry" MR_GROUP_COUNT==1 & MR_DUP_COUNT==0
        }
        scalar MR_GROUP_COUNT=0
        quietly weights_tdnn x, x0matrix(MR_TARGETS) s1(3) s2(5) dim(1)
        test_assert "TDNN weights prepare groups once" MR_GROUP_COUNT==1
        scalar MR_GROUP_COUNT=0
        quietly np_tdnn y, xcols(x) x0matrix(MR_TARGETS) s1(3) s2(5) dim(1)
        test_assert "TDNN estimation reuses weight groups in covariance and posting" MR_GROUP_COUNT==1
    }
    local rc=_rc
    capture adopath - "`MR_SHADOW'"
    quietly do "src/_rdistnn_mata.ado"
    foreach file in group_count _rdistnn_blocks _rdistnn_rank_weights {
        capture erase "`MR_SHADOW'/`file'.ado"
    }
    rmdir "`MR_SHADOW'"
    if `rc' exit `rc'
end
_MR_group_reuse_fixture
program drop _MR_group_reuse_fixture

test_section "reordered saved target profiles reuse cached geometry"
clear
set obs 24
gen double x=mod(_n,7)
gen double y=_n
matrix MR_ORDERED=(0\1\2)
quietly np_dnn y, xcols(x) x0matrix(MR_ORDERED) s(4)
matrix MR_REORDERED=(2\0\2\1)
// A deliberately unavailable distance engine distinguishes cache selection
// from a coincidentally equal recomputed result.
global rdistnn_distance_engine invalid
capture quietly rdistnn_predict, x0matrix(MR_REORDERED) sefit
local MR_REORDER_RC=_rc
macro drop rdistnn_distance_engine
test_assert "reordered and repeated retained targets need no distance evaluation" `MR_REORDER_RC'==0
