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

test_section "scaled norms retain collectively representable subnormal coordinates"
foreach engine in mata plugin {
    global rdistnn_distance_engine `engine'
    foreach metric in l1_scaled l2_scaled {
        mata {
            dr_eta=rdistnn_power2(-1074)
            dr_X=(0,0\dr_eta,dr_eta\2*dr_eta,2*dr_eta\3*dr_eta,3*dr_eta)
            dr_D=rdistnn_distance_matrix(dr_X,J(2,2,0),"`metric'",(2,2))
            dr_want=(0\1\2\3)
            if("`metric'"=="l2_scaled") dr_want=(0\1\1\2)
            st_numscalar("dr_collective",all(dr_D:/dr_eta:==dr_want))
            dr_one=rdistnn_distance_matrix(dr_eta,0,"`metric'",1/.5000000000000001)
            st_numscalar("dr_above_half",dr_one==dr_eta)
        }
        test_assert "`engine'/`metric': analytic final-rounded subnormal norm" scalar(dr_collective)==1
        test_assert "`engine'/`metric': restoration preserves next-above-half mass" scalar(dr_above_half)==1
        test_assert "`engine'/`metric': requested numerical engine actually ran" "$rdistnn_last_distance_engine"=="`engine'"
    }
}
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

// Binary inputs make the final distance exactly 2^1021, although both
// the raw and half-sized positive differences exceed Stata's finite ceiling.
test_section "scaled Mata rows retain quarter-rescued differences"
local dr_quarter_engine "$rdistnn_distance_engine"
global rdistnn_distance_engine mata
mata {
    dr_quarter_x=rdistnn_power2(1022)
    dr_quarter_target=-3*dr_quarter_x
    dr_quarter_want=rdistnn_power2(1021)
    dr_quarter=rdistnn_scaled_difference(dr_quarter_x,dr_quarter_target,.125)
    st_numscalar("dr_quarter_scalar",dr_quarter==dr_quarter_want)
    dr_quarter_l1=rdistnn_distance_matrix(dr_quarter_x,dr_quarter_target,"l1_scaled",8)
    dr_quarter_l2=rdistnn_distance_matrix(dr_quarter_x,dr_quarter_target,"l2_scaled",8)
    st_numscalar("dr_quarter_norms",dr_quarter_l1==dr_quarter_want & dr_quarter_l2==dr_quarter_want)
    st_numscalar("dr_half_restore",rdistnn_distance_power2(.5000000000000001,-1074)==rdistnn_power2(-1074))
}
test_assert "Mata: quarter subtraction restores the finite scaled difference" scalar(dr_quarter_scalar)
test_assert "Mata: scaled L1 and L2 preserve quarter-rescued exact distances" scalar(dr_quarter_norms)
test_assert "Mata: half-split restoration keeps next-above-half subnormal mass" scalar(dr_half_restore)
test_assert "Mata: quarter rescue actually uses the requested engine" "$rdistnn_last_distance_engine"=="mata"
global rdistnn_distance_engine "`dr_quarter_engine'"

// Direct products lie just above or below subnormal half-way points. A
// normalized mantissa product can round away that bit before restoration.
// Zero padding must not turn a one-coordinate norm into another operation.
test_section "scaled norms preserve directly rounded active coordinates"
foreach engine in mata plugin {
    global rdistnn_distance_engine `engine'
    foreach metric in l1_scaled l2_scaled {
        mata {
            dr_eta=rdistnn_power2(-1074)
            dr_inverse=(.16666666666666669,.3)
            dr_count=(3,5)
            dr_direct_ok=dr_scalar_ok=dr_padded_ok=1
            for(dr_case=1;dr_case<=2;dr_case++) {
                dr_raw=dr_count[dr_case]*dr_eta
                dr_inv=dr_inverse[dr_case]
                dr_direct_ok=dr_direct_ok & (dr_raw*dr_inv==dr_eta)
                dr_D=rdistnn_distance_matrix(dr_raw,0,"`metric'",1/dr_inv)
                dr_scalar_ok=dr_scalar_ok & (dr_D==dr_eta)
                dr_D=rdistnn_distance_matrix((0,dr_raw,0),J(2,3,0),"`metric'",(1,1/dr_inv,1))
                dr_padded_ok=dr_padded_ok & all(dr_D:==dr_eta)
            }
            st_numscalar("dr_direct_products",dr_direct_ok)
            st_numscalar("dr_active_scalar",dr_scalar_ok)
            st_numscalar("dr_active_padded",dr_padded_ok)
            // The second product rounds to zero but is genuinely active.
            // Its contribution is below normal accumulator precision, so the
            // complete norm must retain the first rounded coordinate bound.
            dr_D=rdistnn_distance_matrix((3*dr_eta,dr_eta),(0,0),"`metric'",(1/dr_inverse[1],1e300))
            st_numscalar("dr_coordinate_floor",dr_D==dr_eta)
        }
        test_assert "`engine'/`metric': direct rounding fixtures straddle half-way boundaries" scalar(dr_direct_products)==1
        test_assert "`engine'/`metric': one active coordinate keeps its direct rounded product" scalar(dr_active_scalar)==1
        test_assert "`engine'/`metric': padded zeros preserve the one-active-coordinate norm" scalar(dr_active_padded)==1
        test_assert "`engine'/`metric': multiple active coordinates preserve the rounded coordinate bound" scalar(dr_coordinate_floor)==1
        test_assert "`engine'/`metric': direct-rounding route actually ran" "$rdistnn_last_distance_engine"=="`engine'"
    }
}
macro drop rdistnn_distance_engine
