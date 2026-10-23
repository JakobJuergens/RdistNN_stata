*! version 0.1.0  2026-10-22
version 17
// Empty output means known absent capability only. Execution failures propagate.
capture mata: mata drop rdistnn_distance_plugin()
mata:
real matrix rdistnn_distance_plugin(real matrix X, real matrix X0,
    string scalar distance, real matrix geometry)
{
    real scalar rc,i,n,m,api
    real matrix out
    string scalar metric,scale_name
    string rowvector names
    if(st_global("rdistnn_plugin_ok")!="1") return(J(0,0,.))
    names=st_tempname(1)
    st_matrix(names[1],J(1,1,.))
    rc=_stata("plugin call rdistnn_plugin, distance_api "+names[1],1)
    out=st_matrix(names[1])
    st_matrix(names[1],J(0,0,.))
    // Unknown operations are the documented capability-absence protocol.
    if(rc==198) return(J(0,0,.))
    if(rc) rdistnn_fail(rc,"Native distance capability probe failed. Resolve the reported execution/resource error; use rdistnn_distance_engine mata to select Mata explicitly.")
    if(rows(out)!=1 | cols(out)!=1) rdistnn_fail(498,"Native distance capability returned malformed output. Reinstall the current plugin or select Mata explicitly.")
    api=out[1,1]
    if(api==0) return(J(0,0,.))
    if(api!=1) rdistnn_fail(498,"Native distance capability returned an invalid API version. Reinstall the current plugin or select Mata explicitly.")
    n=rows(X); m=rows(X0)
    names=st_tempname(4)
    st_matrix(names[1],X); st_matrix(names[2],X0)
    scale_name="none"
    if(distance=="l1_scaled" | distance=="l2_scaled") {
        st_matrix(names[3],1:/geometry)
        scale_name=names[3]
    }
    metric=(distance=="l1" | distance=="l1_scaled") ? "l1" : "l2"
    st_matrix(names[4],J(n,m,.))
    rc=_stata("plugin call rdistnn_plugin, distance_batch "+metric+" "+
        names[1]+" "+names[2]+" "+scale_name+" "+names[4],1)
    out=st_matrix(names[4])
    for(i=1;i<=4;i++) st_matrix(names[i],J(0,0,.))
    if(rc) rdistnn_fail(rc,"Native distance execution failed. Resolve the reported execution/resource error; select Mata explicitly if needed.")
    if(rows(out)!=n | cols(out)!=m) rdistnn_fail(498,"Native distance execution returned incompatible dimensions. Reinstall the current plugin or select Mata explicitly.")
    if(hasmissing(out) | any(out:<0)) rdistnn_fail(498,"Native distance execution returned invalid distances. Rescale covariates; if this persists, report a reproducible example.")
    return(out)
}
end
