// Analytic extreme L2 cases; forced engine markers prevent fallback masquerading.
version 17
quietly do "src/_rdistnn_mata.ado"
do "tests/test_framework.do"
test_section "L2 range: finite huge/tiny norms and exact neighbor decisions"
foreach engine in mata plugin {
    global rdistnn_distance_engine `engine'
    foreach h in 1e200 1e-200 {
        foreach scaled in 0 1 {
            mata {
                dr_h=strtoreal("`h'")
                dr_X=(3,4\4,3\0,1\0,0\5,0):*dr_h
                dr_T=J(3,2,0)
                dr_G=J(0,0,.)
                dr_metric="l2"
                if(`scaled') {
                    dr_X=dr_X:*(2,4)
                    dr_G=(2,4)
                    dr_metric="l2_scaled"
                }
                dr_D=rdistnn_distance_matrix(dr_X,dr_T,dr_metric,dr_G)
                st_numscalar("dr_error",max(abs(dr_D:/dr_h:-(5\5\1\0\5))))
                st_numscalar("dr_ranks",all(rdistnn_ranks_from_dmat(dr_D):==(3\4\2\1\5)))
                dr_S=rdistnn_distance_matrix(dr_X,dr_T[1,],dr_metric,dr_G)
                st_numscalar("dr_scalar",all(dr_S:==dr_D[,1]))
                dr_sub=rdistnn_distance_matrix(dr_X[(3\1),],dr_T,dr_metric,dr_G)
                st_numscalar("dr_subset",all(dr_sub:==dr_D[(3\1),]))
            }
            test_assert "`engine'/`h'/`scaled': analytic 3-4-5 distances" scalar(dr_error)<2e-14
            test_assert "`engine'/`h'/`scaled': stable neighbor ranks" scalar(dr_ranks)==1
            test_assert "`engine'/`h'/`scaled': scalar equals batch" scalar(dr_scalar)==1
            test_assert "`engine'/`h'/`scaled': frozen subset geometry" scalar(dr_subset)==1
            test_assert "`engine'/`h'/`scaled': actual requested route" "$rdistnn_last_distance_engine"=="`engine'"
        }
    }
    mata {
        dr_X=(1e200,1e-200,0\1e-200,0,0\0,0,0)
        dr_D=rdistnn_distance_matrix(dr_X,J(1,3,0),"l2",J(0,0,.))
        st_numscalar("dr_mixed",all(dr_D:==(1e200\1e-200\0)))
        // Both operands are below Stata's missing-value boundary. Their raw
        // difference is outside its finite range; inverse scaling rescues it.
        dr_D=rdistnn_distance_matrix(8e307,-8e307,"l2_scaled",4)
        st_numscalar("dr_rescue",dr_D/4e307)
    }
    test_assert "`engine': mixed magnitudes and zeros" scalar(dr_mixed)==1
    test_assert "`engine': scaled difference range recovery" abs(scalar(dr_rescue)-1)<2e-15
    capture mata: dr_D=rdistnn_distance_matrix((8e307,8e307),(0,0),"l2",J(0,0,.))
    test_assert "`engine': true Stata result overflow remains explicit" _rc!=0
}
// A later true overflow or Stata missing input must leave the host output intact.
matrix DR_X=(0,0\8e307,8e307)
matrix DR_T=(0,0)
matrix DR_OUT=J(2,1,73)
capture plugin call rdistnn_plugin, distance_batch l2 DR_X DR_T none DR_OUT
local dr_rc=_rc
test_assert "native final overflow rejects batch atomically" `dr_rc'!=0 & DR_OUT[1,1]==73 & DR_OUT[2,1]==73
matrix DR_X[2,1]=.a
capture plugin call rdistnn_plugin, distance_batch l2 DR_X DR_T none DR_OUT
local dr_rc=_rc
test_assert "native Stata missing value rejects batch atomically" `dr_rc'!=0 & DR_OUT[1,1]==73 & DR_OUT[2,1]==73
macro drop rdistnn_distance_engine

// Inputs and final norms are valid Stata doubles, but intermediates cross its
// 2^1023 boundary. Both engines must take identical exceptional arithmetic.
test_section "L2 range: common Stata intermediate boundary"
foreach engine in mata plugin {
    global rdistnn_distance_engine `engine'
    mata {
        dr_X=(7.385580696697091e153,7.827350851331151e153\1.076170170452885e154,0)
        dr_D=rdistnn_distance_matrix(dr_X,J(2,2,0),"l2",J(0,0,.))
        st_numscalar("dr_boundary_values",all(dr_D:==(1.0761701704528852e154\1.076170170452885e154)))
        st_numscalar("dr_boundary_ranks",all(rdistnn_ranks_from_dmat(dr_D):==(2\1)))
        dr_X=(8.000000000000001e307\8e307)
        dr_D=rdistnn_distance_matrix(dr_X,J(2,1,-7.9e307),"l2_scaled",10)
        st_numscalar("dr_boundary_scaled",all(dr_D:==(1.5900000000000001e307\1.59e307)))
        st_numscalar("dr_boundary_scaled_ranks",all(rdistnn_ranks_from_dmat(dr_D):==(2\1)))
    }
    test_assert "`engine': common squared-sum recovery" scalar(dr_boundary_values)==1
    test_assert "`engine': common squared-sum ranks" scalar(dr_boundary_ranks)==1
    test_assert "`engine': half-scaled subtraction" scalar(dr_boundary_scaled)==1
    test_assert "`engine': half-scaled subtraction ranks" scalar(dr_boundary_scaled_ranks)==1
    test_assert "`engine': boundary route actually ran" "$rdistnn_last_distance_engine"=="`engine'"
}
macro drop rdistnn_distance_engine

test_section "L2 range: singleton ordinary rows have no exceptional indices"
foreach engine in mata plugin {
    global rdistnn_distance_engine `engine'
    capture mata: dr_single=rdistnn_distance_matrix((3,4),(0,0\3,4),"l2",J(0,0,.))
    local dr_rc=_rc
    test_assert "`engine': one ordinary row accepts empty recovery index" `dr_rc'==0
    if `dr_rc'==0 {
        mata: st_numscalar("dr_single_ok",all(dr_single:==(5,0)))
        test_assert "`engine': singleton ordinary and zero norms" scalar(dr_single_ok)==1
    }
    capture mata: dr_single=rdistnn_distance_matrix((6,16),(0,0\6,16),"l2_scaled",(2,4))
    local dr_rc=_rc
    test_assert "`engine': one scaled row accepts empty recovery index" `dr_rc'==0
    if `dr_rc'==0 {
        mata: st_numscalar("dr_single_ok",all(dr_single:==(5,0)))
        test_assert "`engine': singleton frozen scales and zero norms" scalar(dr_single_ok)==1
    }
}
macro drop rdistnn_distance_engine
