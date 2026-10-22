version 17
do "tests/test_framework.do"
test_section "balanced native deletion transport"
mata {
    PT_V=(1::17)
    PT_F=rdistnn_native_frame(PT_V)
    PT_FLAT=vec(PT_F')
    st_numscalar("PT_SHAPE",rows(PT_F)==5 & cols(PT_F)==4 & all(PT_FLAT[1..17]:==PT_V) & all(PT_FLAT[18..20]:==0))
    PT_R=((1::6),(6::1)); PT_S=((1::6),sin((1::6))); PT_T=(1\-.5\0\2\-1); PT_G=(1\1\2\2\3\3)
    PT_GROUP_MATA=rdistnn_loo_from_table(PT_R,PT_S,PT_T,0,"mata",PT_G)
    PT_A=rdistnn_native_loo(PT_R,PT_S,PT_T,PT_G)
    PT_B=rdistnn_block_loo_plugin(PT_R,PT_S,PT_T,PT_G)
    st_numscalar("PT_API",rdistnn_native_packed_api()==1)
    st_numscalar("PT_BLOCK",rows(PT_A)==6 & cols(PT_A)==4 & !hasmissing(PT_A) & all(PT_A:==PT_B))
    PT_A=rdistnn_native_loo(PT_R,PT_S,PT_T)
    PT_B=rdistnn_loo_plugin(PT_R,PT_S,PT_T)
    PT_C=rdistnn_loo_prefix(PT_R,PT_S,PT_T)
    st_numscalar("PT_RANK",rows(PT_A)==6 & cols(PT_A)==4 & !hasmissing(PT_A) & all(PT_A:==PT_B) & all(PT_A:==PT_C))
    PT_F=rdistnn_native_frame(vec(PT_S))
    PT_U=rdistnn_native_unframe(PT_F,6,2)
    st_numscalar("PT_RESTORE",all(PT_U:==PT_S))
}
test_assert "balanced frames retain row-major physical order and zero padding" PT_SHAPE
test_assert "fresh native binary advertises packed deletion capability" PT_API
test_assert "grouped bridge executes packed native multi-target and multi-component deletion" PT_BLOCK
test_assert "rank bridge executes packed native multi-target and multi-component deletion" PT_RANK
test_assert "unframing restores logical column-major matrices" PT_RESTORE

// A recognized packed request failure cannot retry a legacy raw operation.
mata: mata drop rdistnn_native_loo()
mata:
real matrix rdistnn_native_loo(real matrix ranks, real matrix signal, real colvector table, | real colvector groups)
{
    return(J(0,0,.))
}
end
capture mata: PT_A=rdistnn_loo_from_table(PT_R,PT_S,PT_T,0,"plugin")
local packed_fail_rc = _rc
test_assert "forced native reports recognized packed failure" `packed_fail_rc' == 498
mata {
    PT_A=rdistnn_loo_from_table(PT_R,PT_S,PT_T,0,"auto")
    st_numscalar("PT_FALLBACK",st_global("rdistnn_last_jackknife_engine")=="mata" & all(PT_A:==PT_C))
}
test_assert "auto uses Mata after recognized packed failure" PT_FALLBACK
capture mata: PT_A=rdistnn_loo_from_table(PT_R,PT_S,PT_T,0,"plugin",PT_G)
local packed_block_fail_rc = _rc
test_assert "forced grouped native reports recognized packed failure" `packed_block_fail_rc' == 498
mata {
    PT_A=rdistnn_loo_from_table(PT_R,PT_S,PT_T,0,"auto",PT_G)
    st_numscalar("PT_BLOCK_FALLBACK",st_global("rdistnn_last_jackknife_engine")=="mata" & all(PT_A:==PT_GROUP_MATA))
}
test_assert "auto grouped path uses Mata after recognized packed failure" PT_BLOCK_FALLBACK
local packed_prior_use_plugin "$rdistnn_use_plugin"
global rdistnn_use_plugin 1
capture mata: PT_A=rdistnn_loo_from_table(PT_R,PT_S,PT_T,0,"auto")
local packed_required_fail_rc = _rc
test_assert "required native auto mode reports recognized packed failure" `packed_required_fail_rc' == 498
global rdistnn_use_plugin "`packed_prior_use_plugin'"
quietly do "src/_rdistnn_native_transport.ado"
