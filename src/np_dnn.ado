*! version 0.0.1  2026-04-05  Jakob Juergens
*! Nonparametric regression at a point via distributional nearest neighbors

program define np_dnn, eclass
    version 16

    if replay() {
        if "`e(cmd)'" != "np_dnn" error 301
        Display `0'
        exit
    }

    syntax varname(numeric) [if] [in] , ///
        x(varlist numeric)               ///  covariates
        x0(numlist)                      ///  evaluation point
        [                                ///
        k(integer 5)                     ///  # nearest neighbors
        s(integer 0)                     ///  kernel order (0 = auto)
        Distance(string)                 ///  l2 (default) | mahalanobis
        Level(cilevel)                   ///
        ]

    marksample touse
    markout `touse' `x'

    local outcome `varlist'
    local xvars   `x'
    local kx : word count `xvars'

    quietly count if `touse'
    if r(N) == 0 { di as error "No observations."; exit 2000 }
    local n = r(N)

    local nx0 : word count `x0'
    if `nx0' != `kx' {
        di as error "x0() must have `kx' elements."
        exit 198
    }

    if "`distance'" == "" local distance "l2"
    if `s' == 0  local s = floor(`n'^((`kx' + 1) / (`kx' + 2)))
    if `s' < `k' {
        di as error "s(`s') < k(`k')."
        exit 198
    }

    quietly run "_rdistnn_mata.ado"

    local x0_mata ""
    foreach v of local x0 {
        local x0_mata `x0_mata' `v',
    }
    local x0_mata = substr("`x0_mata'", 1, length("`x0_mata'") - 1)

    mata {
        _X   = st_data(., "`xvars'", "`touse'")
        _y   = st_data(., "`outcome'", "`touse'")
        _x0  = (`x0_mata')
        _w   = rdistnn_weights(_X, _x0, `k', `s')

        _mu_hat = sum(_w :* _y)

        // LOO jackknife SE
        _n     = rows(_w)
        _denom = 1 :- _w
        _safe  = (_denom :> 1e-10)
        _mu_loo = J(_n, 1, _mu_hat)
        _mu_loo[selectindex(_safe)] = ///
            (_mu_hat :- _w[selectindex(_safe)] :* _y[selectindex(_safe)]) ///
            :/ _denom[selectindex(_safe)]

        _se_jk = sqrt((_n - 1) / _n * sum((_mu_loo :- _mu_hat):^2))
        _neff  = 1 / sum(_w:^2)

        st_numscalar("r(_mu_hat)", _mu_hat)
        st_numscalar("r(_se_jk)",  _se_jk)
        st_numscalar("r(_neff)",   _neff)
    }

    tempname b V
    matrix `b' = (r(_mu_hat))
    matrix `V' = (r(_se_jk)^2)
    matrix colnames `b' = mu
    matrix rownames `V' = mu
    matrix colnames `V' = mu

    ereturn post `b' `V', obs(`n') esample(`touse')
    ereturn scalar neff   = r(_neff)
    ereturn scalar k      = `k'
    ereturn scalar s      = `s'
    ereturn scalar kx     = `kx'
    ereturn local  x0     "`x0'"
    ereturn local  xvars  "`xvars'"
    ereturn local  cmd    "np_dnn"
    ereturn local  cmdline "np_dnn `0'"

    Display, level(`level')
end

program define Display
    syntax [, Level(cilevel)]
    di ""
    di as text "DNN Nonparametric Regression"
    di as text "Evaluation point x0: " as result "`e(x0)'"
    di as text "k = " as result e(k) as text "  s = " as result e(s) ///
       as text "  neff = " as result %6.2f e(neff)
    di ""
    ereturn display, level(`level')
end
