version 17
do "tests/test_framework.do"
test_section "original-signal block deletion numerical regressions"
foreach engine in naive mata plugin {
    mata {
        BN_R=(1::3); BN_X=(0\0\1); BN_T=(1\0)
        BN_L=rdistnn_block_loo(BN_R,(2^54\1\0),BN_T,0,"`engine'",BN_X)
        st_numscalar("BN_SURVIVOR",BN_L[1]==1 & BN_L[2]==2^54 & BN_L[3]==2^53)
        st_numscalar("BN_ENGINE",st_global("rdistnn_last_jackknife_engine")=="`engine'")
    }
    test_assert "`engine': dominant duplicate deletion retains exact finite survivor" BN_SURVIVOR
    test_assert "`engine': block deletion executes the requested engine" BN_ENGINE
    mata {
        BN_L=rdistnn_block_loo((1::4),(2^54\1\2^53\2^53),(0\1\-1),0,"`engine'",(0\0\1\2))
        BN_V=rdistnn_covariance(BN_L,3/4,1)
        st_numscalar("BN_VARIATION",all(BN_L:==(0\0\.5\.5)) & BN_V[1,1]==3/16)
        BN_ETA=rdistnn_power2(-1074)
        BN_L=rdistnn_block_loo(BN_R,(BN_ETA\0\0),BN_T,0,"`engine'",BN_X)
        st_numscalar("BN_SUBNORMAL",BN_ETA>0 & all(BN_L:==(0\BN_ETA\0)))
        BN_L=rdistnn_block_loo(BN_R,(BN_ETA\0\0),(rdistnn_power2(1000)\0),0,"`engine'",BN_X)
        st_numscalar("BN_AMPLIFIED",all(BN_L:==(0\rdistnn_power2(-74)\rdistnn_power2(-75))))
        BN_L=rdistnn_loo_from_table(BN_R,(BN_ETA\0\0),(.5\.5),0,"`engine'")
        st_numscalar("BN_RANK_SUBNORMAL",all(BN_L:==0))
    }
    test_assert "`engine': signed zero-mass tables preserve exact nonzero covariance" BN_VARIATION
    test_assert "`engine': genuine subnormal signals retain rounded-product deletion semantics" BN_SUBNORMAL
    test_assert "`engine': unnormalized mass recovers amplified subnormal contributions" BN_AMPLIFIED
    test_assert "`engine': rank-only and block paths share exceptional numerical origins" BN_RANK_SUBNORMAL
    mata {
        BN_R=(1\9\3\7\5\2\4\6\8\10)
        BN_X=(J(5,1,0)\(1::5))
        BN_T=(rdistnn_power2(112)\0\rdistnn_power2(56)\0\1\0\-rdistnn_power2(56)\0\-rdistnn_power2(112))
        BN_Z=(0\1\1\1\1\J(5,1,0))
        BN_L=rdistnn_block_loo(BN_R,BN_Z,BN_T,0,"`engine'",BN_X)
        st_numscalar("BN_OBSERVATION_ORDER",BN_L[10]==4/5)
        st_numscalar("BN_OBSERVATION_ENGINE",st_global("rdistnn_last_jackknife_engine")=="`engine'")
    }
    test_assert "`engine': interleaved blocks retain the public observation-order mass" BN_OBSERVATION_ORDER
    test_assert "`engine': interleaved mass regression executes the requested engine" BN_OBSERVATION_ENGINE
    mata {
        BN_R=(1\5\2\4\3\6)
        BN_X=(J(5,1,0)\1)
        BN_T=(rdistnn_power2(112)\rdistnn_power2(56)\1\-rdistnn_power2(56)\-rdistnn_power2(112))
        BN_Z=(0\1\1\1\1\0)
        BN_L=rdistnn_block_loo(BN_R,BN_Z,BN_T,0,"`engine'",BN_X)
        st_numscalar("BN_MEMBER_ORDER",BN_L[6]==4/5)
        st_numscalar("BN_MEMBER_ENGINE",st_global("rdistnn_last_jackknife_engine")=="`engine'")
    }
    test_assert "`engine': contiguous supplied blocks retain original member-order masses" BN_MEMBER_ORDER
    test_assert "`engine': member-order regression executes the requested engine" BN_MEMBER_ENGINE
}
