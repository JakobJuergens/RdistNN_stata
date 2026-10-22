version 17
do "tests/test_framework.do"
do "tests/helpers/finite-matrix.do"
test_section "duplicate raw covariate block averaging and actual deletions"
mata: st_numscalar("BLOCK_COMPARE_MATCH",rdistnn_test_finite_close((0\1),(0\1),0))
mata: st_numscalar("BLOCK_COMPARE_PARTIAL",rdistnn_test_finite_close((0\.),(0\1),0))
mata: st_numscalar("BLOCK_COMPARE_ALL",rdistnn_test_finite_close((.\.),(0\1),0))
scalar BLOCK_COMPARE_SHAPE=1
capture mata: st_numscalar("BLOCK_COMPARE_SHAPE",rdistnn_test_finite_close((0,0),(0\0),0))
local compare_shape_rc=_rc
test_assert "finite comparator accepts matching represented values" BLOCK_COMPARE_MATCH
test_assert "finite comparator rejects a partially missing actual result" BLOCK_COMPARE_PARTIAL==0
test_assert "finite comparator rejects an entirely missing actual result" BLOCK_COMPARE_ALL==0
test_assert "finite comparator rejects wrong shape without aborting" `compare_shape_rc'==0 & BLOCK_COMPARE_SHAPE==0
// Independent reference: enumerate equal raw rows after explicitly removing
// each observation, and average the selected rank table over each surviving
// group. Neither production block helpers nor mean-signal corrections are used.
capture mata: mata drop rdistnn_test_block_weights()
capture mata: mata drop rdistnn_test_block_deletions()
mata:
real matrix rdistnn_test_block_weights(real matrix X,real matrix R,real colvector table)
{
    real scalar i,j,t,n
    real colvector members
    real matrix W
    n=rows(X); W=J(n,cols(R),.)
    for(i=1;i<=n;i++) {
        members=J(n,1,0)
        for(j=1;j<=n;j++) members[j]=all(X[i,]:==X[j,])
        for(t=1;t<=cols(R);t++) W[i,t]=quadsum(table[select(R[,t],members)])/sum(members)
    }
    return(W)
}
real matrix rdistnn_test_block_deletions(real matrix X,real matrix R,
    real matrix S,real colvector table)
{
    real scalar i,t,n,q
    real colvector kept
    real matrix RR,W,L
    n=rows(X); q=cols(S); L=J(n,cols(R)*q,.)
    for(i=1;i<=n;i++) {
        kept=selectindex((1::n):!=i)
        RR=R[kept,]:-(R[kept,]:>J(n-1,1,1)*R[i,])
        W=rdistnn_test_block_weights(X[kept,],RR,table)
        for(t=1;t<=cols(R);t++) L[i,((t-1)*q+1)..(t*q)]=(S[kept,]'*W[,t])'
    }
    return(L)
}
BX=(0,0\0,0\0,0\1,0\1,0\-1,0\0,1\2,1)
BS=(1,4\7,2\-2,5\4,-1\8,3\3,6\5,8\-1,2)
BT=(0,0\.5,0\0,0\1,1)
BP=(5\3\8\1\7\4\2\6)
BERR=1; BLERR=1; BCERR=1; BPERM=1; BENGINE=1
end
foreach metric in l1 l2 l1_scaled l2_scaled mahalanobis {
    foreach method in exact asymp {
        foreach k in 1 2 {
            mata {
                BG=rdistnn_fit_geometry(BX,"`metric'","sd",J(1,0,.),1e-12)
                BD=rdistnn_distance_matrix(BX,BT,"`metric'",BG)
                BR=rdistnn_localization_ranks(BD,BX)
                BW=rdistnn_weights_from_rmat(BR,`k',5,"`method'",BX)
                BF=rdistnn_weights_from_ranks((1::8),`k',5,"`method'")
                BREF=rdistnn_test_block_weights(BX,BR,BF)
                BERR=BERR & rdistnn_test_finite_close(BW,BREF,2e-15)
                BDEL=rdistnn_weights_from_ranks((1::7),`k',5,"`method'")
                BLREF=rdistnn_test_block_deletions(BX,BR,BS,BDEL)
                BVREF=rdistnn_covariance(BLREF,7/8,1)
                BPR=BR[BP,]; BPW=rdistnn_weights_from_rmat(BPR,`k',5,"`method'",BX[BP,])
                BPERM=BPERM & rdistnn_test_finite_close(BPW,BW[BP,],2e-13)
            }
            foreach engine in naive mata plugin {
                mata {
                    BL=rdistnn_loo_from_table(BR,BS,BDEL,0,"`engine'",BX)
                    BLERR=BLERR & rdistnn_test_finite_close(BL,BLREF,2e-13)
                    BV=rdistnn_covariance(rdistnn_loo_from_table(BR,BS,BDEL,1,"`engine'",BX),7/8,1)
                    BCERR=BCERR & rdistnn_test_finite_close(BV,BVREF,5e-12)
                    BPL=rdistnn_loo_from_table(BPR,BS[BP,],BDEL,0,"`engine'",BX[BP,])
                    BPERM=BPERM & rdistnn_test_finite_close(BPL,BL[BP,],2e-13)
                    BENGINE=BENGINE & st_global("rdistnn_last_jackknife_engine")=="`engine'"
                }
            }
        }
    }
}
mata: st_numscalar("BLOCK_W",BERR); st_numscalar("BLOCK_L",BLERR); st_numscalar("BLOCK_V",BCERR); st_numscalar("BLOCK_P",BPERM); st_numscalar("BLOCK_ENGINE",BENGINE)
test_assert "raw-X block weights equal independent averaging across metrics, exact/asymp, k=1/2" BLOCK_W
test_assert "naive/Mata/native deletion arrays equal independently regrouped delete-one reference" BLOCK_L
test_assert "all components and joint target covariance equal explicit deletion covariance" BLOCK_V
test_assert "row permutation preserves weights and aligned deletion arrays" BLOCK_P
test_assert "every forced deletion engine actually executed, including native plugin" BLOCK_ENGINE & "$rdistnn_plugin_ok"=="1"

// Explicit grouping metadata also respects arbitrary supplied ranks, where
// equal raw rows need not occupy contiguous positions. Rank-only behavior is
// otherwise unchanged.
mata {
    IX=(0\1\0\2\0\1\2\3); IR=(1::8)
    IS=(1,4\7,2\-2,5\4,-1\8,3\3,6\5,8\-1,2)
    IT=(.8\-.2\.3\.1\0\0\0)
    IL=rdistnn_test_block_deletions(IX,IR,IS,IT)
    IW=(.7\.2\-.1\.3\-.1\0\0\0)
    IWR=rdistnn_test_block_weights(IX,IR,IW)
    st_numscalar("BLOCK_INTERW",rdistnn_test_finite_close(rdistnn_block_weights(IR,IW,IX),IWR,1e-15))
}
foreach engine in naive mata plugin {
    mata: st_numscalar("BLOCK_INTER",rdistnn_test_finite_close(rdistnn_loo_from_table(IR,IS,IT,0,"`engine'",IX),IL,1e-13))
    test_assert "`engine': explicitly grouped interleaved ranks retain actual deletion shifts" BLOCK_INTER
}
test_assert "explicitly grouped interleaved rank weights preserve group mass" BLOCK_INTERW

findfile "_rdistnn_tdnn_mata.ado"
quietly run "`r(fn)'"
mata: BTERR=1; BTLERR=1; BTVERR=1
foreach method in exact asymp {
    foreach normalize in 0 1 {
        mata {
            BF=rdistnn_tdnn_table(8,2,5,2,"`method'",`normalize')
            BW=rdistnn_tdnn_weights(BR,2,5,2,"`method'",`normalize',BX)
            BTERR=BTERR & rdistnn_test_finite_close(BW,rdistnn_test_block_weights(BX,BR,BF),2e-15)
            BDEL=rdistnn_tdnn_table(7,2,5,2,"`method'",`normalize')
            BLREF=rdistnn_test_block_deletions(BX,BR,BS,BDEL)
        }
        foreach engine in naive mata plugin {
            global rdistnn_jackknife_engine `engine'
            mata {
                BL=rdistnn_loo_from_table(BR,BS,BDEL,0,"`engine'",BX)
                BTLERR=BTLERR & rdistnn_test_finite_close(BL,BLREF,2e-13)
                BV=rdistnn_tdnn_cov(BR,BS,2,5,2,"`method'",`normalize',BX)
                BTVERR=BTVERR & rdistnn_test_finite_close(BV,rdistnn_covariance(BLREF,7/8,1),5e-12)
            }
        }
    }
}
mata: st_numscalar("BLOCK_TW",BTERR); st_numscalar("BLOCK_TL",BTLERR); st_numscalar("BLOCK_TV",BTVERR)
test_assert "signed and unnormalized TDNN block weights retain reference block mass" BLOCK_TW
test_assert "all TDNN deletion engines implement actual-observation removal" BLOCK_TL
test_assert "signed TDNN full joint covariance agrees with explicit deletion" BLOCK_TV
macro drop rdistnn_jackknife_engine

// A finite deletion may require an unrepresentable unweighted difference.
// Its ordinary variance overflows, so test raw deletion arrays separately.
foreach engine in naive mata plugin {
    mata {
        RX=J(5,1,0); RY=(8e307\8e307\8e307\8e307\-8e307)
        RL=rdistnn_loo_from_table((1::5),RY,rdistnn_exact_table(4,1,2),0,"`engine'",RX)
        st_numscalar("BLOCK_CORRECTION_RANGE",rdistnn_test_finite_close(RL:/(J(4,1,4e307)\8e307),J(5,1,1),2e-15))
    }
    test_assert "`engine': finite huge deletion survives overflowing unweighted correction" BLOCK_CORRECTION_RANGE
    mata {
        RX=(0\0\0\1\2); RY=(8e307\8e307\8e307\-1e307\0)
        RL=rdistnn_loo_from_table((1::5),RY,rdistnn_exact_table(4,1,2),0,"`engine'",RX)
        st_numscalar("BLOCK_MERGE_RANGE",!hasmissing(RL) & rdistnn_test_finite_close(RL[1..3,1]:/6.5e307,J(3,1,1),2e-15))
    }
    test_assert "`engine': zero block correction merges without intermediate overflow" BLOCK_MERGE_RANGE
}

// Every identical X row must have equal weight; jackknife reduces to the usual
// sample-mean deletions even if exact support cuts through the block.
mata {
    BA=J(8,2,0); R=1::8; S=(1::8)
    W=rdistnn_weights_from_rmat(R,1,7,"exact",BA)
    L=rdistnn_loo_from_table(R,S,rdistnn_exact_table(7,1,7),0,"plugin",BA)
    st_numscalar("BLOCK_ALL",rdistnn_test_finite_close(W,J(8,1,1/8),1e-15) & rdistnn_test_finite_close(L,(36:-S)/7,1e-14))
    st_numscalar("BLOCK_ALLV",rdistnn_test_finite_close(rdistnn_jk_cov_multi_from_ranks(R,S,1,7,"exact",BA),variance(S)/8,1e-14))
    // The support boundary enters the entire 3-row first block. The remaining
    // extreme group is unused and must not affect origin selection or zeros.
    XS=(0\0\0\1\1\2\3\4)
    SS=(1\5\-2\8e307\8e307\-8e307\8e307\-8e307)
    RS=(1::8); TS=rdistnn_exact_table(7,1,7)
    LS=rdistnn_loo_from_table(RS,SS,TS,1,"plugin",XS)
    LF=rdistnn_test_block_deletions(XS,RS,SS,TS)
    st_numscalar("BLOCK_SUPPORT",!hasmissing(LS) & !hasmissing(LF) & rdistnn_test_finite_close(rdistnn_covariance(LS,7/8,1),rdistnn_covariance(LF,7/8,1),1e-14))
    Z=J(8,1,2^100)
    st_numscalar("BLOCK_ZERO",rdistnn_test_finite_close(rdistnn_jk_cov_multi_from_ranks(RS,Z,1,7,"exact",XS),0,0))
    BLOCK_ETA=rdistnn_power2(-1074)
    TW=rdistnn_block_weights((1::3),(BLOCK_ETA\BLOCK_ETA\0),(0\0\1))
    st_numscalar("BLOCK_TINY",BLOCK_ETA>0 & rdistnn_test_finite_close(TW,(BLOCK_ETA\BLOCK_ETA\0),0))
    HW=rdistnn_block_weights((1::3),(8e307\8e307\0),(0\0\1))
    st_numscalar("BLOCK_HUGE",rdistnn_test_finite_close(HW,(8e307\8e307\0),0))
}
test_assert "all-identical X produces equal weights and usual leave-one means" BLOCK_ALL
test_assert "all-identical X produces usual sample-mean variance" BLOCK_ALLV
test_assert "support expansion includes full duplicate block and excludes extreme unused groups" BLOCK_SUPPORT
test_assert "constant supported signal preserves exact zero-variance status" BLOCK_ZERO
test_assert "block averaging retains nonzero subnormal mass and exact unused zeros" BLOCK_TINY
test_assert "block averaging recovers finite means from overflowing intermediate mass" BLOCK_HUGE
global rdistnn_plugin_ok 0
capture mata: rdistnn_loo_from_table((1::8),J(8,1,0),rdistnn_exact_table(7,1,7),1,"plugin",J(8,1,0))
local zero_native_rc=_rc
global rdistnn_plugin_ok 1
test_assert "forced native zero-signal deletion cannot bypass an unavailable plugin" `zero_native_rc'==498

clear
set obs 8
mata: st_store(.,st_addvar("double",("bx1","bx2","by","bd")),(BX,BS))
gen long bid=_n
gen double bone=1
gen double bzero=0
matrix BX0=(0,0\.5,0\0,0)
foreach family in np cate clate tdnn plqr {
    local cmd "`family'_dnn"
    local xopt "xcols(bx1 bx2)"
    local extra "s(5)"
    local input "by"
    local targets "x0matrix(BX0)"
    if "`family'"=="cate" local extra "`extra' pseudooutcome"
    if "`family'"=="clate" local extra "`extra' pseudooutcome denominator(bd)"
    if "`family'"=="tdnn" {
        local cmd "np_tdnn"
        local extra "s1(2) s2(5)"
    }
    if "`family'"=="plqr" {
        local xopt "x(bx1 bx2)"
        local extra "s(5) d(by bd) md(bzero bzero) g(bzero) thetaeval(0 0)"
        local input "bone"
        local targets "x0(0 0)"
    }
    global rdistnn_jackknife_engine plugin
    quietly `cmd' `input', `xopt' `targets' `extra'
    matrix BB=e(estimates)
    matrix BV=e(estimate_vcov)
    test_assert "`family': public fit retained new tie policy and used native deletions" "`e(tie_rule)'"=="distance_then_covariates_with_duplicate_block_averaging" & e(model_version)==2 & "$rdistnn_last_jackknife_engine"=="plugin"
    quietly rdistnn_results, type(weights)
    matrix BWW=r(data)
    mata: st_numscalar("BLOCK_EXTRACT",rdistnn_test_finite_close(st_matrix("BWW")[1..3,5],J(3,1,st_matrix("BWW")[1,5]),1e-15))
    test_assert "`family': retained weight extraction agrees inside duplicate block" BLOCK_EXTRACT
    if "`family'"=="clate" matrix BCV=e(component_vcov)
    preserve
    gsort -bid
    global rdistnn_jackknife_engine mata
    quietly `cmd' `input', `xopt' `targets' `extra'
    mata: st_numscalar("BLOCK_FIT",rdistnn_test_finite_close(st_matrix("e(estimates)"),st_matrix("BB"),1e-12)); st_numscalar("BLOCK_COV",rdistnn_test_finite_close(st_matrix("e(estimate_vcov)"),st_matrix("BV"),1e-11))
    test_assert "`family': public estimates are permutation invariant across Mata/native" BLOCK_FIT
    test_assert "`family': public full covariance is permutation invariant across Mata/native" BLOCK_COV
    if "`family'"=="clate" {
        mata: st_numscalar("BLOCK_COMP",rdistnn_test_finite_close(st_matrix("e(component_vcov)"),st_matrix("BCV"),1e-11))
        test_assert "CLATE paired component covariance is permutation invariant" BLOCK_COMP
    }
    restore
    tempfile snapshot
    quietly estimates save "`snapshot'", replace
    preserve
    clear
    mata: mata clear
    quietly estimates use "`snapshot'"
    if "`family'"=="plqr" quietly rdistnn_predict, x0(0 0) sefit post
    else quietly rdistnn_predict, x0matrix(BX0) sefit post
    quietly do "tests/helpers/finite-matrix.do"
    mata: st_numscalar("BLOCK_SAVED",rdistnn_test_finite_close(st_matrix("e(estimates)"),st_matrix("BB"),1e-12) & rdistnn_test_finite_close(st_matrix("e(estimate_vcov)"),st_matrix("BV"),1e-11))
    test_assert "`family': cold saved-model predictions retain block estimates and covariance" BLOCK_SAVED
    restore
}
macro drop rdistnn_jackknife_engine
quietly np_dnn by, xcols(bx1 bx2) x0matrix(BX0) tuningmode(grid) sgrid(4 5) sdefault(5)
matrix BTABLE=e(tuning_diagnostics)
estimates store duplicate_grid
quietly rdistnn_candidate, candidate(1) post
quietly rdistnn_results, type(weights)
matrix BWW=r(data)
mata: st_numscalar("BLOCK_GRID",rdistnn_test_finite_close(rdistnn_grid_weight_diagnostics((st_matrix("BWW")[1..8,5],st_matrix("BWW")[9..16,5],st_matrix("BWW")[17..24,5])),st_matrix("BTABLE")[1,],1e-12))
test_assert "candidate diagnostics summarize actual target weights" BLOCK_GRID
quietly estimates restore duplicate_grid
quietly np_dnn_plot, type(weights) view(profile) allcandidates targets(2 1) nodraw replace
matrix BPD=r(plot_data)
test_assert "candidate weight profiles retain target identifiers" colsof(BPD)==7 & rowsof(BPD)==32
quietly np_dnn_plot, type(weights) view(cumulative) targets(2 1) nodraw replace
matrix BPD=r(plot_data)
test_assert "cumulative profiles retain the requested target-specific weights" colsof(BPD)==7 & rowsof(BPD)==16 & BPD[1,7]==2 & BPD[9,7]==1
quietly np_tdnn by, xcols(bx1 bx2) x0matrix(BX0) s1(2) s2(5)
quietly np_tdnn_plot, type(weights) targets(2 1) nodraw replace
matrix BPD=r(plot_data)
test_assert "signed TDNN profiles reconstruct target-specific averaged weights" colsof(BPD)==3 & rowsof(BPD)==16
// Old fit schemas may retain old estimates. They cannot reconstruct weights or
// covariance under the new policy until refitted.
capture program drop rdistnn_test_old_model
program define rdistnn_test_old_model, eclass
    ereturn scalar model_version=1
end
rdistnn_test_old_model
capture quietly rdistnn_predict, x0(0 0) sefit
local oldrc=_rc
test_assert "old saved schemas require refitting before prediction" `oldrc'==498
