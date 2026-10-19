*! version 0.1.0  2026-10-18
version 17
// One matrix bridge for scalar/batch L1/L2, including fitted inverse scales.
capture mata: mata drop rdistnn_distance_plugin()
mata:
real matrix rdistnn_distance_plugin(real matrix X, real matrix X0,
    string scalar distance, real matrix geometry)
{
    real scalar rc,i,api,n,m
    real matrix out
    string scalar cached,metric,scale_name
    string rowvector names
    if(st_global("rdistnn_plugin_ok")!="1") return(J(0,0,.))
    cached=st_global("rdistnn_distance_plugin_api")
    if(cached!="0" & cached!="1") {
        names=st_tempname(1)
        st_matrix(names[1],J(1,1,.))
        rc=_stata("plugin call rdistnn_plugin, distance_api "+names[1],1)
        api=0
        if(rc==0) {
            out=st_matrix(names[1])
            if(rows(out)==1 & cols(out)==1) api=(out[1,1]==1)
        }
        st_matrix(names[1],J(0,0,.))
        st_global("rdistnn_distance_plugin_api",api ? "1" : "0")
        cached=st_global("rdistnn_distance_plugin_api")
    }
    if(cached!="1") return(J(0,0,.))
    n=rows(X); m=rows(X0)
    names=st_tempname(4)
    st_matrix(names[1],X)
    st_matrix(names[2],X0)
    scale_name="none"
    if(distance=="l1_scaled" | distance=="l2_scaled") {
        st_matrix(names[3],1:/geometry)
        scale_name=names[3]
    }
    metric=(distance=="l1" | distance=="l1_scaled") ? "l1" : "l2"
    st_matrix(names[4],J(n,m,.))
    rc=_stata("plugin call rdistnn_plugin, distance_batch "+metric+" "+
        names[1]+" "+names[2]+" "+scale_name+" "+names[4],1)
    out=J(0,0,.)
    if(rc==0) {
        out=st_matrix(names[4])
        if(rows(out)!=n | cols(out)!=m) out=J(0,0,.)
        else if(hasmissing(out) | any(out:<0)) out=J(0,0,.)
    }
    for(i=1;i<=4;i++) st_matrix(names[i],J(0,0,.))
    return(out)
}
end
