version 17
do "tests/test_framework.do"
quietly do "src/_rdistnn_mata.ado"
test_section "CATE population signal contracts"
// The normalized Robinson signal identifies the oracle effect but has linear
// propensity drift away from e0=.5. AIPW retains the true effect with true arms.
mata:
_fn_y=(0\0\0\1); _fn_d=_fn_y
_fn_e=.26
_fn_robin=mean((_fn_d:-_fn_e):*(_fn_y:-.25)/(_fn_e*(1-_fn_e)))
_fn_aipw=mean(J(4,1,1)+_fn_d:*(_fn_y:-1)/_fn_e-(1:-_fn_d):*_fn_y/(1-_fn_e))
st_numscalar("_fn_score",abs(_fn_robin-.1875/(_fn_e*(1-_fn_e)))<1e-14 & _fn_aipw==1)
end
test_assert "Robinson and AIPW retain their distinct population contracts" scalar(_fn_score)
