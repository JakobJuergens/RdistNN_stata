version 17
do "tests/test_framework.do"
test_section "weighted underflow numerical origins"
mata {
    WO_ETA=rdistnn_power2(-1074)
    WO_R=(1::3); WO_X=(0\0\1)
    WO_S=(1.5,1.5\0,1.5\0,1.5)
    WO_T=(WO_ETA\WO_ETA)
    WO_EXPECT=(0,4\2,4\2,4):*WO_ETA
}
foreach engine in naive mata plugin {
    foreach center in 0 1 {
        mata {
            WO_L=rdistnn_loo_from_table(WO_R,WO_S,WO_T,`center',"`engine'")
            st_numscalar("WO_RANK",!hasmissing(WO_L) & all(WO_L:==WO_EXPECT))
            st_numscalar("WO_ENGINE",st_global("rdistnn_last_jackknife_engine")=="`engine'")
            WO_L=rdistnn_loo_from_table(WO_R,WO_S,WO_T,`center',"`engine'",WO_X)
            st_numscalar("WO_BLOCK",!hasmissing(WO_L) & all(WO_L:==WO_EXPECT))
            WO_W=1.25*rdistnn_power2(-1022)
            WO_L=rdistnn_loo_from_table(WO_R,J(3,1,1.5+epsilon(1)),(WO_W\-(WO_W-WO_ETA)),`center',"`engine'")
            st_numscalar("WO_SIGNED",!hasmissing(WO_L) & all(WO_L:==WO_ETA))
        }
        test_assert "`engine' center=`center': zero survivors and constants retain rounded products" WO_RANK
        test_assert "`engine' center=`center': requested actual engine executes" WO_ENGINE
        test_assert "`engine' center=`center': duplicate paths share weighted origins" WO_BLOCK
        test_assert "`engine' center=`center': nonzero subnormal signed mass preserves constants" WO_SIGNED
    }
    foreach units in 1 4 {
        mata {
            WO_UNITS=`units'
            WO_L=rdistnn_loo_from_table((1::4),J(4,1,rdistnn_power2(1000)),(1\-1\(WO_UNITS*WO_ETA)),0,"`engine'",J(4,1,0))
            WO_VALUE=(WO_UNITS==4 ? 3*rdistnn_power2(-74) : 0)
            st_numscalar("WO_AVERAGE",!hasmissing(WO_L) & all(WO_L:==WO_VALUE))
        }
        test_assert "`engine': amplified rounded duplicate coefficient with mass `units' eta" WO_AVERAGE
    }
}
mata {
    st_numscalar("WO_TAIL",rdistnn_signal_level((1.5\0\0),(1\WO_ETA))==1.5)
    st_numscalar("WO_MASS_ZERO",rdistnn_signal_level(J(3,1,1.5),(1\-1))==1.5)
}
test_assert "ordinary tiny tails preserve their observed origin" WO_TAIL
test_assert "exact zero-mass normal signed tables retain centering" WO_MASS_ZERO
