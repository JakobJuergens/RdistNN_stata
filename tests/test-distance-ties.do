version 17
do "tests/test_framework.do"
do "tests/helpers/finite-matrix.do"
test_section "exact represented distance block averaging"
clear
input double x y
-1 1
1 9
-2 4
2 8
3 5
end
quietly weights_dnn x, x0(0) s(4)
matrix DISTANCE_W=r(weights)
// The first two distinct covariates are equidistant. Their combined mass
// is one; their individual ranks have different raw table weights.
mata: st_numscalar("DISTANCE_EQUAL",rdistnn_test_finite_close(st_matrix("DISTANCE_W"),(0.5\0.5\0\0\0),2e-15))
test_assert "different covariates at equal metric distance share total rank mass" DISTANCE_EQUAL
quietly np_dnn y, xcols(x) x0(0) s(4)
scalar DISTANCE_MEAN=e(b)[1,1]
test_assert "estimator uses the averaged actual distance weights" abs(DISTANCE_MEAN-5)<1e-13

do "tests/helpers/distance-reference.do"
clear
input double x1 x2 y den
1 0 2 1.1
-1 0 8 1.5
0 1 -3 1.7
0 -1 4 .8
2 0 7 1.3
-2 0 1 2
0 2 -1 .7
3 1 5 1.6
end
gen long distance_id=_n
gen double distance_one=1
gen double distance_zero=0
matrix DISTANCE_TARGETS=(0,0\1,0\-1,0)
mata: DT_X=st_data(.,"x1 x2"); DT_S=st_data(.,"y den"); DT_TARGETS=st_matrix("DISTANCE_TARGETS")
foreach metric in l1 l2 l1_scaled l2_scaled mahalanobis {
    foreach distance_engine in mata plugin {
        global rdistnn_distance_engine `distance_engine'
        local scales ""
        if inlist("`metric'","l1_scaled","l2_scaled") local scales "scalerule(custom) scalevec(1 2)"
        quietly weights_dnn x1 x2, x0matrix(DISTANCE_TARGETS) k(2) s(6) distance(`metric') `scales'
        matrix DT_R=r(ranks)
        matrix DT_D=r(distances)
        matrix DT_RETURNED_GROUPS=r(groups)
        matrix DT_W=r(weights)
        mata {
            DT_D=st_matrix("DT_D"); DT_R=st_matrix("DT_R"); DT_W=st_matrix("DT_W")
            DT_RETG=st_matrix("DT_RETURNED_GROUPS")
            DT_GROUP_OK=1
            for(DT_GI=1;DT_GI<=8;DT_GI++) for(DT_GJ=1;DT_GJ<=8;DT_GJ++) {
                DT_GROUP_OK=DT_GROUP_OK & all((DT_RETG[DT_GI,]:==DT_RETG[DT_GJ,]):==(DT_D[DT_GI,]:==DT_D[DT_GJ,]))
            }
            st_numscalar("DT_RETURNED_GROUPS",DT_GROUP_OK)
            DT_REF=rdistnn_test_subset_weights(DT_D,2,6)
            DT_LREF=rdistnn_test_subset_deletions(DT_D,DT_S,2,6)
            DT_VREF=(7/8)*cross(DT_LREF:-mean(DT_LREF),DT_LREF:-mean(DT_LREF))
            st_numscalar("DT_OK",rdistnn_test_finite_close(DT_W,DT_REF,3e-15))
        }
        test_assert "`distance_engine'/`metric': actual weights match independent subset selection with tied boundaries" DT_OK
        test_assert "`distance_engine'/`metric': returned groups equal exact represented-distance partition" DT_RETURNED_GROUPS
        test_assert "`distance_engine'/`metric': requested distance engine executed" "$rdistnn_last_distance_engine"=="`distance_engine'"
        foreach deletion_engine in naive mata plugin {
            mata {
                DT_L=rdistnn_loo_from_ranks(DT_R,DT_S,2,6,"exact",0,
                    rdistnn_weights_from_ranks((1::7),2,6,"exact"),DT_D)
            }
            // Explicit engine argument exercises the common dispatcher, rather
            // than allowing a coincidental default to satisfy the assertion.
            mata {
                DT_L=rdistnn_loo_from_table(DT_R,DT_S,rdistnn_exact_table(7,2,6),0,"`deletion_engine'",DT_D)
                DT_C=rdistnn_loo_from_table(DT_R,DT_S,rdistnn_exact_table(7,2,6),1,"`deletion_engine'",DT_D)
                st_numscalar("DT_L_OK",rdistnn_test_finite_close(DT_L,DT_LREF,3e-13))
                st_numscalar("DT_V_OK",rdistnn_test_finite_close(rdistnn_covariance(DT_C,7/8,1),DT_VREF,8e-12))
            }
            test_assert "`distance_engine'/`metric'/`deletion_engine': actual-observation deletions match subset enumeration" DT_L_OK
            test_assert "`distance_engine'/`metric'/`deletion_engine': centered paired and cross-target covariance matches enumeration" DT_V_OK
            test_assert "`distance_engine'/`metric'/`deletion_engine': selected deletion engine ran" "$rdistnn_last_jackknife_engine"=="`deletion_engine'"
        }
    }
}
macro drop rdistnn_distance_engine

// Signed, nonunit mass tables keep each exact block's signed mass, including
// a block crossing the nominal support boundary and target-specific groups.
mata {
    DT_D=rdistnn_distance_matrix(DT_X,DT_TARGETS,"l1")
    DT_R=rdistnn_localization_ranks(DT_D)
    DT_TABLE=(.7\-.3\.2\0\-.1\0\0)
    DT_LREF=rdistnn_test_distance_deletions(DT_D,DT_R,DT_S,DT_TABLE)
}
foreach engine in naive mata plugin {
    mata {
        DT_L=rdistnn_loo_from_table(DT_R,DT_S,DT_TABLE,0,"`engine'",DT_D)
        st_numscalar("DT_SIGNED",rdistnn_test_finite_close(DT_L,DT_LREF,1e-13))
        DT_EXPLICIT=rdistnn_distance_loo((1::8),DT_S,DT_TABLE,0,"`engine'",(0\1\0\2\0\1\2\3))
        DT_EREF=rdistnn_test_distance_deletions((0\1\0\2\0\1\2\3),(1::8),DT_S,DT_TABLE)
        st_numscalar("DT_EXPLICIT",rdistnn_test_finite_close(DT_EXPLICIT,DT_EREF,1e-13))
    }
    test_assert "`engine': signed unnormalized tables retain exact-distance block mass after actual deletion" DT_SIGNED
    test_assert "`engine': authoritative interleaved supplied ranks exercise direct group kernel" DT_EXPLICIT
    test_assert "`engine': interleaved group calculation ran selected engine" "$rdistnn_last_jackknife_engine"=="`engine'"
}

foreach engine in naive mata plugin {
    global rdistnn_jackknife_engine `engine'
    quietly weights_tdnn x1 x2, x0matrix(DISTANCE_TARGETS) s1(3) s2(5) dim(2) distance(l1)
    matrix DT_TW=r(weights)
    mata {
        DT_TREF=-1.5*rdistnn_test_subset_weights(DT_D,1,3)+2.5*rdistnn_test_subset_weights(DT_D,1,5)
        DT_TLREF=-1.5*rdistnn_test_subset_deletions(DT_D,DT_S[,1],1,3)+2.5*rdistnn_test_subset_deletions(DT_D,DT_S[,1],1,5)
        st_numscalar("DT_TW_OK",rdistnn_test_finite_close(st_matrix("DT_TW"),DT_TREF,3e-15))
        DT_TVREF=(7/8)*cross(DT_TLREF:-mean(DT_TLREF),DT_TLREF:-mean(DT_TLREF))
    }
    test_assert "`engine': signed TDNN weights equal independent combination of subset learners" DT_TW_OK
    quietly np_tdnn y, xcols(x1 x2) x0matrix(DISTANCE_TARGETS) s1(3) s2(5) dim(2) distance(l1)
    mata: st_numscalar("DT_TV_OK",rdistnn_test_finite_close(st_matrix("e(estimate_vcov)"),DT_TVREF,8e-12))
    test_assert "`engine': signed TDNN covariance equals independent subset deletions" DT_TV_OK
    test_assert "`engine': TDNN estimator used requested deletion engine" "$rdistnn_last_jackknife_engine"=="`engine'"
}
macro drop rdistnn_jackknife_engine

quietly weights_dnn x1 x2, x0matrix(DISTANCE_TARGETS) s(8) distance(l1)
matrix DT_BOUNDARY=r(weights)
test_assert "nominal one-rank support expands to all four equally nearest covariates" DT_BOUNDARY[1,1]==.25 & DT_BOUNDARY[2,1]==.25 & DT_BOUNDARY[3,1]==.25 & DT_BOUNDARY[4,1]==.25 & DT_BOUNDARY[5,1]==0
mata {
    DT_NEAR=(1\1+2^-48\2\3)
    DT_NR=rdistnn_localization_ranks(DT_NEAR)
    DT_NW=rdistnn_weights_from_rmat(DT_NR,1,4,"exact",DT_NEAR)
    st_numscalar("DT_NEAR",all(DT_NW:==(1\0\0\0)))
    DT_PURE=rdistnn_weights_from_rmat((1::4),1,4,"exact")
    st_numscalar("DT_PURE",all(DT_PURE:==(1\0\0\0)))
}
test_assert "unequal represented distances remain distinct despite being extremely close" DT_NEAR
test_assert "rank-only primitives retain strict supplied-rank semantics without accompanying distances" DT_PURE
mata: st_numscalar("DT_GROUP_INTER",all(rdistnn_distance_groups((1::4),(1\2\1\2)):==(1\2\1\2)))
test_assert "distance groups use exact equality independently of supplied rank adjacency" DT_GROUP_INTER
foreach engine in naive mata plugin {
    mata: st_numscalar("DT_STRICT_LOO",all(rdistnn_loo_from_table((1::4),(1\9\4\8),(1\0\0),0,"`engine'"):==(9\1\1\1)))
    test_assert "`engine': rank-only calls retain strict supplied-rank deletions" DT_STRICT_LOO
}

quietly weights_dnn x1 x2, x0matrix(DISTANCE_TARGETS) s(6) distance(l2)
matrix DT_ORIGINAL=r(weights)
gen double reflected_x1=-x2
gen double reflected_x2=x1
matrix DISTANCE_REFLECTED=(0,0\0,1\0,-1)
quietly weights_dnn reflected_x1 reflected_x2, x0matrix(DISTANCE_REFLECTED) s(6) distance(l2)
mata: st_numscalar("DT_REFLECT",rdistnn_test_finite_close(st_matrix("r(weights)"),st_matrix("DT_ORIGINAL"),2e-15))
test_assert "signed coordinate reflection and permutation preserve distance weights" DT_REFLECT

// Saved models retain target-specific distance groups across cold reuse.
quietly np_dnn y, xcols(x1 x2) x0matrix(DISTANCE_TARGETS) tuningmode(grid) sgrid(4 6) sdefault(6) distance(l1)
matrix DT_FIT=e(estimates)
matrix DT_V=e(estimate_vcov)
matrix DT_GROUPS=e(model_groups)
test_assert "new fits retain schema3 and exact distance averaging policy" e(model_version)==3 & "`e(tie_rule)'"=="exact_distance_block_averaging_v1" & rowsof(DT_GROUPS)==8 & colsof(DT_GROUPS)==3
test_assert "different targets retain different exact distance groups" DT_GROUPS[1,1]==DT_GROUPS[2,1] & DT_GROUPS[1,2]!=DT_GROUPS[2,2]
tempfile distance_cold
quietly estimates save "`distance_cold'", replace
preserve
clear
mata: mata clear
quietly estimates use "`distance_cold'"
quietly rdistnn_predict, x0matrix(DISTANCE_TARGETS) sefit post
do "tests/helpers/finite-matrix.do"
mata: st_numscalar("DT_COLD",rdistnn_test_finite_close(st_matrix("e(estimates)"),st_matrix("DT_FIT"),1e-13) & rdistnn_test_finite_close(st_matrix("e(estimate_vcov)"),st_matrix("DT_V"),3e-12) & all(st_matrix("e(model_groups)"):==st_matrix("DT_GROUPS")))
test_assert "cold prediction preserves fitted distances groups estimates and covariance" DT_COLD
quietly rdistnn_predict, x0(0 1) sefit post
matrix DT_NEW_GROUPS=e(model_groups)
test_assert "new target prediction caches its own distance groups" colsof(DT_NEW_GROUPS)==1 & rowsof(DT_NEW_GROUPS)==8
restore
