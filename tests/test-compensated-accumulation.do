version 17
quietly do "src/_rdistnn_mata.ado"
do "tests/test_framework.do"
test_section "compensated weighted sums and retained deletion corrections"
mata:
_ca_A=2^54
_ca_R=(1::4)
_ca_S=((_ca_A\1\0\-_ca_A),J(4,1,7))
_ca_T=J(3,1,1)
_ca_expected=(-_ca_A,21\0,21\1,21\_ca_A,21)
end
foreach engine in mata naive plugin {
    capture noisily mata: _ca_L=rdistnn_loo_from_table(_ca_R,_ca_S,_ca_T,0,"`engine'")
    local rc=_rc
    test_assert "`engine': deletion engine actually available" `rc'==0
    if `rc'==0 {
        mata: st_numscalar("_ca_equal",all(_ca_L:==_ca_expected))
        test_assert "`engine': preserve small contribution across final merge" scalar(_ca_equal)
    }
}
capture mata: st_numscalar("_ca_sum",rdistnn_sum((_ca_A\1\-_ca_A)))
local rc=_rc
test_assert "shared compensated sum helper available" `rc'==0
if `rc'==0 {
    test_assert "signed reduction retains unit contribution" scalar(_ca_sum)==1
    mata: st_numscalar("_ca_large",rdistnn_sum((1\2^200\1\-2^200)))
    test_assert "large incoming term retains both unit contributions" scalar(_ca_large)==2
    mata: _ca_X=((_ca_A\1\-_ca_A),J(3,1,1)); _ca_W=(J(3,1,1/3),(0\1\0))
    mata: st_numscalar("_ca_weighted",all(rdistnn_weighted_signals(_ca_X,_ca_W):==(1/3,1\1,1)))
    test_assert "paired weighted signals use the shared arithmetic" scalar(_ca_weighted)
}
foreach engine in mata plugin naive {
    mata: _ca_L=rdistnn_loo_from_table((1::2),(2\3),1,0,"`engine'")
    mata: st_numscalar("_ca_edge",all(_ca_L:==(3\2)))
    test_assert "`engine': two-observation deletions" scalar(_ca_edge)
    mata: _ca_L=rdistnn_loo_from_table((1::4),(1\2\3\2^800),(1\0\0),0,"`engine'")
    mata: st_numscalar("_ca_edge",all(_ca_L:==(2\1\1\1)))
    test_assert "`engine': unused huge values cannot contaminate support" scalar(_ca_edge)
}
mata:
_ca_n=4101
_ca_z=J(_ca_n,1,0)
_ca_z[(1\2048\2050\4101)]=(2^54\1\-2^54\1)
_ca_R=((1::_ca_n),(_ca_n::1))
_ca_signal=(_ca_z,J(_ca_n,1,7))
_ca_table=J(_ca_n-1,1,1)
_ca_native=rdistnn_loo_from_table(_ca_R,_ca_signal,_ca_table,0,"plugin")
_ca_mata=rdistnn_loo_from_table(_ca_R,_ca_signal,_ca_table,0,"mata")
st_numscalar("_ca_blocks",all(_ca_native:==_ca_mata))
st_numscalar("_ca_anchor",all(_ca_mata[2049,(1,3)]:==2))
end
test_assert "large paired deletions agree exactly across block boundaries" scalar(_ca_blocks)
test_assert "block-boundary deletions retain both unit contributions" scalar(_ca_anchor)

// Native quad products must retain the same separately rounded-product contract.
mata: st_numscalar("_ca_product",rdistnn_weighted_signals((1+2^-27\1),(1-2^-27\-1)))
test_assert "quad weighted signals use separately rounded products" scalar(_ca_product)==0
mata: _ca_missing=rdistnn_weighted_signals((.\1),(1\1))
mata: st_numscalar("_ca_missing1",missing(_ca_missing))
test_assert "weighted signals propagate a missing input" scalar(_ca_missing1)
mata: _ca_missing=rdistnn_weighted_signals((. ,2\1,3),(1,1\1,1))
mata: st_numscalar("_ca_missing2",all(_ca_missing[1,]:>=.) & all(_ca_missing[2,]:==5))
test_assert "missing signal component cannot omit rows from valid components" scalar(_ca_missing2)
mata: _ca_overflow=rdistnn_weighted_signals((2^800\1),(2^800\1))
mata: st_numscalar("_ca_overflowed",missing(_ca_overflow))
test_assert "weighted product overflow propagates missing" scalar(_ca_overflowed)

clear
set obs 4
gen double _ca_x=_n
gen double _ca_y=cond(_n==2,2^54,cond(_n==4,-2^54,cond(_n==3,1,0)))
quietly np_tdnn _ca_y, xcols(_ca_x) x0(0) s1(1) s2(4) nose
matrix _ca_b=e(b)
quietly weights_tdnn _ca_x, x0(0) s1(1) s2(4)
mata: st_numscalar("_ca_tdnn",all(st_matrix("_ca_b"):==rdistnn_weighted_signals(st_data(.,"_ca_y"),st_matrix("r(weights)"))))
test_assert "TDNN posted point uses shared signed-weight arithmetic" scalar(_ca_tdnn)
