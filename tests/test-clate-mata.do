// Paired CLATE evaluator: literal R fixture and independent subset deletions.
version 17
quietly do "src/_rdistnn_mata.ado"
capture confirm file "src/_rdistnn_clate_mata.ado"
if !_rc quietly do "src/_rdistnn_clate_mata.ado"
do "tests/test_framework.do"

test_section "CLATE Mata: four-observation analytic fixture"
mata {
    _cl_R = (1::4)
    _cl_S = (2,1 \ 1,-1 \ 4,2 \ 3,1)
    _cl_C = _cl_V = _cl_CV = _cl_status = J(0,0,.)
    _cl_t = J(1,0,.)
    rdistnn_clate_fit(_cl_R, _cl_S, 1, 2, "exact", 1, 1e-8,
        _cl_C, _cl_t, _cl_V, _cl_CV, _cl_status)
    st_numscalar("_cl_comp", max(abs(_cl_C - (2,.5))))
    st_numscalar("_cl_ratio", abs(_cl_t[1]-4))
    st_numscalar("_cl_var", abs(_cl_V[1,1]-34))
    st_numscalar("_cl_covar", max(abs(_cl_CV - (.5,.5 \ .5,.75))))
    st_numscalar("_cl_statusok", max(abs(_cl_status - (0,0,0,0,1))))
}
test_assert "components are weighted means, never individual ratios" scalar(_cl_comp) < 1e-12
test_assert "ratio is four" scalar(_cl_ratio) < 1e-12
test_assert "linearized paired ratio variance is 34" scalar(_cl_var) < 1e-10
test_assert "paired component covariance retains cross term" scalar(_cl_covar) < 1e-12
test_assert "regular target has Wald inference" scalar(_cl_statusok) == 0

// Literal subset enumeration is independent of DNN weight/deletion helpers.
// For s=3, enumerate all triples and average the first k by frozen rank.
capture mata: mata drop clate_test_subsets()
mata:
real rowvector clate_test_subsets(real colvector ranks, real matrix sig, real scalar k)
{
    real scalar a,b,c,n,count
    real colvector ix, ord
    real rowvector value
    n=rows(ranks); count=0; value=(0,0)
    for(a=1;a<=n-2;a++) for(b=a+1;b<=n-1;b++) for(c=b+1;c<=n;c++) {
        ix=(a\b\c); ord=order(ranks[ix],1)
        value=value+mean(sig[ix[ord[1..k]],]); count++
    }
    return(value/count)
}
end

test_section "CLATE Mata: frozen-rank batch covariance"
mata {
    _cl_R = (1,6 \ 3,4 \ 5,2 \ 6,1 \ 4,3 \ 2,5)
    _cl_S = (1,1 \ 3,2 \ -2,1 \ 4,3 \ 2,2 \ 5,1)
    _cl_idx=(1::6)
    _cl_err=0; _cl_cross=0
    for(_cl_k=1;_cl_k<=2;_cl_k++) {
        rdistnn_clate_fit(_cl_R,_cl_S,_cl_k,3,"exact",1,1e-8,
            _cl_C,_cl_t,_cl_V,_cl_CV,_cl_status)
        _cl_loo=J(6,4,.)
        _cl_full=J(2,2,.)
        for(_cl_j=1;_cl_j<=2;_cl_j++) {
            _cl_full[_cl_j,]=clate_test_subsets(_cl_R[,_cl_j],_cl_S,_cl_k)
            for(_cl_i=1;_cl_i<=6;_cl_i++) {
                _cl_keep=(_cl_idx:!=_cl_i)
                _cl_loo[_cl_i,(2*_cl_j-1)..(2*_cl_j)] = clate_test_subsets(
                    select(_cl_R[,_cl_j],_cl_keep),select(_cl_S,_cl_keep),_cl_k)
            }
        }
        _cl_center=_cl_loo:-mean(_cl_loo)
        _cl_wantCV=5/6*cross(_cl_center,_cl_center)
        _cl_gradient=J(4,2,0)
        for(_cl_j=1;_cl_j<=2;_cl_j++) {
            _cl_gradient[(2*_cl_j-1)..(2*_cl_j),_cl_j] =
                (1 \ -_cl_t[_cl_j])/_cl_full[_cl_j,2]
        }
        _cl_wantV=_cl_gradient'*_cl_wantCV*_cl_gradient
        _cl_err=max((_cl_err,max(abs(_cl_full-_cl_C)),
            max(abs(_cl_wantCV-_cl_CV)),max(abs(_cl_wantV-_cl_V))))
        _cl_cross=max((_cl_cross,abs(_cl_V[1,2])))
    }
    st_numscalar("_cl_subseterr",_cl_err)
    st_numscalar("_cl_cross",_cl_cross)
}
test_assert "general-k paired deletion covariance matches subset enumeration" scalar(_cl_subseterr) < 1e-10
test_assert "batch retains nonzero cross-target covariance" scalar(_cl_cross) > 1e-6

mata {
    _cl_R=(1::4); _cl_S=(2,1 \ 1,-1 \ 4,2 \ 3,1)
    rdistnn_clate_fit(_cl_R,_cl_S,1,2,"asymp",1,1e-8,
        _cl_C,_cl_t,_cl_V,_cl_CV,_cl_status)
    st_numscalar("_cl_asympcomp",max(abs(_cl_C-(17/8,9/16))))
    st_numscalar("_cl_asympvar",abs(_cl_V[1,1]-3075200/177147))
    rdistnn_clate_fit(_cl_R,_cl_S,1,2,"exact",0,1e-8,
        _cl_C,_cl_t,_cl_V,_cl_CV,_cl_status)
    st_numscalar("_cl_nose",all(_cl_V:==.) & all(_cl_CV:==.) &
        _cl_status[1,4]==1 & _cl_status[1,5]==0 & abs(_cl_t[1]-4)<1e-12)
    _cl_S[,1]=4*_cl_S[,2]
    rdistnn_clate_fit(_cl_R,_cl_S,1,2,"exact",1,1e-8,
        _cl_C,_cl_t,_cl_V,_cl_CV,_cl_status)
    st_numscalar("_cl_zerovar",_cl_V[1,1]==0 & _cl_status[1,4]==2 & _cl_status[1,5]==0)
    _cl_S=(2,0 \ 1,0 \ 4,-1 \ 3,2)
    _cl_R=(1,3,4 \ 2,4,3 \ 3,1,2 \ 4,2,1)
    rdistnn_clate_fit(_cl_R,_cl_S,1,3,"exact",1,1e-8,
        _cl_C,_cl_t,_cl_V,_cl_CV,_cl_status)
    st_numscalar("_cl_exception",_cl_t[1]==. & _cl_status[1,1]==1 &
        _cl_status[2,2]==1 & _cl_t[2]<0 & _cl_t[3]<. &
        all(_cl_V[1,]:==.) & !hasmissing(_cl_CV))
    _cl_S=(2,1e-12 \ 1,-1e-12 \ 4,0 \ 3,1)
    rdistnn_clate_fit((1::4),_cl_S,1,3,"exact",0,1e-8,
        _cl_C,_cl_t,_cl_V,_cl_CV,_cl_status)
    st_numscalar("_cl_near",_cl_status[1,3]==1 & _cl_t[1]<. & _cl_t[1]>1e10)
}
test_assert "asymptotic fixture components" scalar(_cl_asympcomp) < 1e-12
test_assert "asymptotic fixture covariance" scalar(_cl_asympvar) < 1e-10
test_assert "nose returns all missing covariance with status" scalar(_cl_nose) == 1
test_assert "zero paired score variance suppresses Wald availability" scalar(_cl_zerovar) == 1
test_assert "zero and signed denominators preserve full target tables" scalar(_cl_exception) == 1
test_assert "near-zero flag never trims the ratio" scalar(_cl_near) == 1

test_section "CLATE Mata: general-k asymptotic analytic fixture"
mata {
    _cl_R=(1::4); _cl_S=(2,1 \ 1,-1 \ 4,2 \ 3,1)
    // Integrated rank weights (47,41,29,11)/128; deleted weights
    // (13,10,4)/27. Expectations were derived with rational arithmetic.
    rdistnn_clate_fit(_cl_R,_cl_S,2,3,"asymp",1,1e-8,
        _cl_C,_cl_t,_cl_V,_cl_CV,_cl_status)
    st_numscalar("_cl_k2comp",max(abs(_cl_C-(71/32,75/128))))
    st_numscalar("_cl_k2cv",max(abs(_cl_CV-(2219/3888,37/72 \ 37/72,571/972))))
    st_numscalar("_cl_k2v",abs(_cl_V[1,1]-114263428096/7688671875))
}
test_assert "general-k asymptotic components match integrated weights" scalar(_cl_k2comp) < 1e-12
test_assert "general-k asymptotic component covariance matches literal deletions" scalar(_cl_k2cv) < 1e-12
test_assert "general-k asymptotic ratio covariance matches rational fixture" scalar(_cl_k2v) < 1e-10

test_section "CLATE Mata: numerical failure statuses"
mata {
    _cl_R=(1::4); _cl_S=(2,1 \ 1,-1 \ 4,2 \ 3,1)
    _cl_scaleerr=0
    for(_cl_scaleindex=1;_cl_scaleindex<=3;_cl_scaleindex++) {
        _cl_scale=(-1e6,1e-6,1e6)[_cl_scaleindex]
        rdistnn_clate_fit(_cl_R,_cl_scale*_cl_S,1,2,"exact",1,1e-8,
            _cl_C,_cl_t,_cl_V,_cl_CV,_cl_status)
        _cl_scaleerr=max((_cl_scaleerr,abs(_cl_t[1]-4),abs(_cl_V[1,1]-34)))
    }
    st_numscalar("_cl_scaleerr",_cl_scaleerr)
    _cl_S=(2e150,1e-150 \ 1e150,-1e-150 \ 4e150,2e-150 \ 3e150,1e-150)
    rdistnn_clate_fit(_cl_R,_cl_S,1,2,"exact",1,1e-8,
        _cl_C,_cl_t,_cl_V,_cl_CV,_cl_status)
    st_numscalar("_cl_overflow",_cl_t[1]<. & _cl_V[1,1]==. &
        _cl_status[1,4]==3 & _cl_status[1,5]==0)
    _cl_S=(2e-200,1 \ 1e-200,1 \ 4e-200,1 \ 3e-200,1)
    rdistnn_clate_fit(_cl_R,_cl_S,1,2,"exact",1,1e-8,
        _cl_C,_cl_t,_cl_V,_cl_CV,_cl_status)
    st_numscalar("_cl_underflow",_cl_status[1,4]==3 &
        _cl_V[1,1]==. & all(_cl_CV:==.))
    _cl_S=(2e150,1e-200 \ 1e150,-1e-200 \ 4e150,2e-200 \ 3e150,1e-200)
    rdistnn_clate_fit(_cl_R,_cl_S,1,2,"exact",0,1e-8,
        _cl_C,_cl_t,_cl_V,_cl_CV,_cl_status)
    st_numscalar("_cl_ratiooverflow",_cl_t[1]==. & _cl_status[1,1]==2)
}
test_assert "common signal scaling preserves ratio and variance" scalar(_cl_scaleerr) < 1e-10
test_assert "overflowing covariance reports numerical failure" scalar(_cl_overflow) == 1
test_assert "underflowing varying components never masquerade as zero variance" scalar(_cl_underflow) == 1
test_assert "overflowing ratio has numerical failure ratio status" scalar(_cl_ratiooverflow) == 1

test_section "CLATE Mata: instrument and paired-score diagnostics"
mata {
    _cl_R=(1::6); _cl_S=(1,1 \ 3,2 \ -2,1 \ 4,3 \ 2,2 \ 5,1)
    _cl_W=rdistnn_weights_from_rmat(_cl_R,1,1,"exact")
    _cl_diag=J(1,25,99)
    _cl_A=(0\1\0\1\0\1); _cl_e=J(6,1,.4)
    _cl_result=rdistnn_clate_diagnostics(_cl_diag,_cl_W,_cl_S,(2),_cl_A,_cl_e)
    st_numscalar("_cl_arms",max(abs(_cl_result[1,12..23]-(1,.5,.5,.5,1,.4,.4,.4,99,1,3,3))))
    // |Y-2D|=(1,1,4,2,2,3); largest five sum 12 of total 13.
    st_numscalar("_cl_contrib",abs(_cl_result[1,25]-12/13))
    _cl_result=rdistnn_clate_diagnostics(_cl_diag,_cl_W,_cl_S,(2),J(0,1,.),J(0,1,.))
    st_numscalar("_cl_noarms",all(_cl_result[1,(12,16,21)]:==0) &
        all(_cl_result[1,(13..15,17..19,22..23)]:==.) & _cl_result[1,20]==99)
    _cl_result=rdistnn_clate_diagnostics(_cl_diag,_cl_W,_cl_S,(.),_cl_A,_cl_e)
    st_numscalar("_cl_undefineddiag",_cl_result[1,24]==0 & _cl_result[1,25]==.)
}
test_assert "overlap uses instrument arm masses ESS and propensity" scalar(_cl_arms) < 1e-12
test_assert "concentration uses paired CLATE score" scalar(_cl_contrib) < 1e-12
test_assert "pseudo route clears unavailable overlap and preserves radius" scalar(_cl_noarms) == 1
test_assert "undefined ratio has unavailable score diagnostic" scalar(_cl_undefineddiag) == 1

// Reload the extension twice without unloading the base library.
quietly do "src/_rdistnn_clate_mata.ado"
quietly do "src/_rdistnn_clate_mata.ado"
mata: st_numscalar("_cl_reload",rows(rdistnn_weights_exact((1::4),1,2)))
test_assert "extension reload preserves base Mata functions" scalar(_cl_reload) == 4
