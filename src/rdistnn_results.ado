*! version 0.1.0  2026-10-18
program define rdistnn_results, rclass
    version 17
    syntax [, TYPE(string) TARGETS(numlist integer >=1) Level(cilevel) SAVING(string asis) REPLACE]
    if !inlist("`e(cmd)'","np_dnn","np_tdnn","cate_dnn","clate_dnn","plqr_dnn") error 301
    if "`type'"=="" local type estimates
    if !inlist("`type'","estimates","vcov","components","status","diagnostics","weights","influence") error 198
    local m=e(m)
    if missing(`m') | `m'<1 | `m'!=floor(`m') {
        di as error "Current fitted results require a positive integer e(m); refit `e(cmd)' with the current package."
        exit 498
    }
    if "`targets'"=="" numlist "1/`m'"
    if "`targets'"=="" local targets "`r(numlist)'"
    foreach t of local targets {
        if `t'>`m' error 198
    }
    local targets : list uniq targets
    tempname table summary fit se V data index coordinates selection
    quietly _rdistnn_report, level(`level')
    matrix `table'=r(table)
    matrix `summary'=r(summary)
    local q=1
    if "`e(cmd)'"=="plqr_dnn" local q=e(kd)
    mata: st_matrix("`selection'",vec(J(`q',1,1)*((strtoreal(tokens("`targets'")):-1):*`q')+(1::`q')*J(1,cols(tokens("`targets'")),1))')
    mata: st_matrix("`table'",st_matrix("`table'")[,st_matrix("`selection'")]); st_matrix("`summary'",st_matrix("`summary'")[st_matrix("`selection'"),])
    matrix `fit'=`table'[1,1...]
    matrix `se'=`table'[2,1...]
    local fields : e(matrices)
    foreach field in estimate_vcov targets {
        if !`: list field in fields' {
            di as error "Current fitted results require e(`field'); refit `e(cmd)' with the current package."
            exit 498
        }
    }
    mata: st_matrix("`V'",st_matrix("e(estimate_vcov)")[st_matrix("`selection'"),st_matrix("`selection'")])
    mata: st_matrix("`fit'",colshape(st_matrix("`fit'"),`q')'); st_matrix("`se'",colshape(st_matrix("`se'"),`q')')
    matrix `coordinates'=e(targets)
    mata: st_matrix("`index'",(floor((st_matrix("`selection'")':-1):/`q'):+1,mod(st_matrix("`selection'")':-1,`q'):+1))
    matrix colnames `index'=target component
    if "`type'"=="estimates" {
        mata: st_matrix("`data'",(st_matrix("`index'"),st_matrix("`coordinates'")[st_matrix("`index'")[,1],],st_matrix("`summary'")))
        local cols "target component"
        forvalues j=1/`=colsof(`coordinates')' {
            local cols "`cols' x`j'"
        }
        matrix colnames `data'=`cols' estimate std_error z pvalue ci_lower ci_upper
    }
    if "`type'"=="vcov" {
        mata: st_matrix("`data'",(st_matrix("`index'")[ceil((1::(rows(st_matrix("`V'"))^2)):/rows(st_matrix("`V'"))),],st_matrix("`index'")[mod((1::(rows(st_matrix("`V'"))^2)):-1,rows(st_matrix("`V'"))):+1,],vec(st_matrix("`V'")') ))
        matrix colnames `data'=row_target row_component column_target column_component covariance
    }
    if inlist("`type'","components","status","diagnostics") {
        local field "`type'"
        if !`:list field in fields' {
            di as error "No `type' are stored for this fit."
            exit 498
        }
        mata: st_matrix("`data'",(strtoreal(tokens("`targets'"))',st_matrix("e(`field')")[strtoreal(tokens("`targets'")),]))
        matrix colnames `data'=target `:colnames e(`field')'
    }
    if inlist("`type'","weights","influence") {
        if e(model_version)!=1 {
            di as error "Observation extraction requires a current saved model; refit `e(cmd)' with the current package."
            exit 498
        }
        if "`type'"=="influence" & !inlist("`e(cmd)'","cate_dnn","clate_dnn") error 198
        foreach file in _rdistnn_mata.ado _rdistnn_clate_mata.ado _rdistnn_model_mata.ado {
            findfile "`file'"
            quietly run "`r(fn)'"
        }
        if "`type'"=="weights" {
            mata: st_matrix("`data'",rdistnn_model_weights(strtoreal(tokens("`targets'"))))
            matrix colnames `data'=observation_id target rank distance weight
        }
        else {
            mata: st_matrix("`data'",rdistnn_model_plot_data("influence",strtoreal(tokens("`targets'"))))
            if "`e(cmd)'"=="clate_dnn" matrix colnames `data'=target observation_id rank distance weight score contribution phi_y phi_d
            else matrix colnames `data'=target observation_id rank distance weight score contribution
        }
    }
    quietly _rdistnn_diagnostics_save, matrix(`data') saving(`saving') `replace'
    if "`e(cmd)'"=="plqr_dnn" {
        tempname theta
        matrix `theta'=e(theta_eval)
        return matrix theta_eval=`theta'
        return local quantity "moment"
    }
    else return local quantity "estimate"
    return matrix fit=`fit'
    return matrix se_fit=`se'
    return matrix vcov=`V'
    return matrix table=`table'
    return matrix summary=`summary'
    return matrix index=`index'
    return matrix data=`data'
    return scalar level=`level'
    return local type "`type'"
end
