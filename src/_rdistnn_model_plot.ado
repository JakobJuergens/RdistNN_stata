*! version 0.1.0  2026-10-18
program define _rdistnn_model_plot, rclass
    version 17
    syntax, TYPE(string) NAME(name) [TARGETS(numlist integer >=1) NLABEL(integer 3) REPLACE NODRAW]
    if "`targets'"=="" {
        if e(m)<=8 numlist "1/`=e(m)'"
        else numlist "1 `=ceil(e(m)/2)' `=e(m)'"
        local targets "`r(numlist)'"
    }
    local ntargets : word count `targets'
    foreach file in _rdistnn_mata.ado _rdistnn_clate_mata.ado _rdistnn_model_mata.ado {
        findfile "`file'"
        quietly run "`r(fn)'"
    }
    tempname plot
    noisily mata: st_matrix("`plot'",rdistnn_model_plot_data("`type'",strtoreal(tokens("`targets'"))))
    local nameopt "name(`name')"
    if "`replace'"!="" local nameopt "name(`name', replace)"
    if "`type'"=="kernel" matrix colnames `plot'=target_id rank distance weight cum_share
    else if "`e(cmd)'"=="clate_dnn" matrix colnames `plot'=target_id obs_id rank distance weight score contribution phi_y phi_d
    else matrix colnames `plot'=target_id obs_id rank distance weight score contribution
    preserve
    clear
    quietly svmat double `plot', names(col)
    if "`type'"=="kernel" {
        twoway (connected weight distance, sort msize(vsmall)), ///
            by(target_id, note("") title("DNN effective kernel")) ///
            xtitle("Distance") ytitle("DNN weight") `nameopt' `nodraw'
    }
    else {
        quietly gen double abs_contribution=abs(contribution)
        quietly gsort target_id -abs_contribution
        quietly by target_id: gen str12 obs_label=string(obs_id,"%9.0g") if _n<=`nlabel' & `nlabel'>0
        quietly sort target_id distance
        twoway (scatter contribution distance, mcolor(black%65) msize(small) mlabel(obs_label) mlabsize(vsmall)), ///
            by(target_id, note("") title("DNN frozen-score contributions")) ///
            yline(0, lpattern(dash) lcolor(gs8)) xtitle("Distance") ///
            ytitle("Weighted centered score") `nameopt' `nodraw'
    }
    restore
    return matrix plot_data=`plot'
    return scalar n_plotted=`ntargets'
    return local view "`type'"
end
