*! version 0.1.5  2026-10-25
version 17
// Balanced private matrix frames avoid expensive default stripes on tall
// native inputs and outputs. All native deletions use this current transport.
foreach fn in frame unframe packed_api loo {
    capture mata: mata drop rdistnn_native_`fn'()
}
mata:
// Physical cells hold a row-major sequence; logical arrays use vec() order.
real matrix rdistnn_native_frame(real colvector payload)
{
    real scalar length,r,c
    length=rows(payload)
    if(length<1 | length>9007199254740991) rdistnn_fail(198,"Native transport requires a nonempty payload with a representable cell count.")
    r=ceil(sqrt(length)); c=ceil(length/r)
    return(rowshape(payload\J(r*c-length,1,0),r))
}

real matrix rdistnn_native_unframe(real matrix frame, real scalar n, real scalar columns)
{
    real colvector payload
    real scalar length
    length=n*columns
    if(n<1 | columns<1 | n!=floor(n) | columns!=floor(columns) |
        length>9007199254740991 | rows(frame)<cols(frame) |
        rows(frame)-cols(frame)>1 | rows(frame)*cols(frame)<length |
        rows(frame)*cols(frame)-length>=rows(frame)) return(J(0,0,.))
    payload=vec(frame')
    if(hasmissing(payload)) return(J(0,0,.))
    if(rows(payload)>length) {
        if(any(payload[(length+1)..rows(payload)]:!=0)) return(J(0,0,.))
    }
    return(rowshape(payload[1..length],columns)')
}

// Probe every call: a program may be rebound to a different binary in-session.
// _stata(...,1) suppresses unsupported-operation output while retaining its rc.
real scalar rdistnn_native_packed_api()
{
    string scalar name
    real scalar rc,api
    real matrix result
    if(st_global("rdistnn_plugin_ok")!="1") return(0)
    name=st_tempname()
    st_matrix(name,J(1,1,.))
    rc=_stata("plugin call rdistnn_plugin, loo_packed_api "+name,1)
    if(!rc) result=st_matrix(name)
    st_matrix(name,J(0,0,.))
    if(rc==198) return(0)
    if(rc) rdistnn_fail(rc,"Native deletion capability probe failed. Resolve the execution/resource error or select rdistnn_jackknife_engine mata explicitly.")
    if(rows(result)!=1 | cols(result)!=1) rdistnn_fail(498,"Native deletion capability returned malformed output. Reinstall the current plugin or select Mata explicitly.")
    api=result[1,1]
    if(api==0) return(0)
    if(api!=1) rdistnn_fail(498,"Native deletion capability returned an invalid version. Reinstall the current plugin or select Mata explicitly.")
    return(1)
}

// Header: version,n,m,q,grouped; then vec(ranks), vec(signal), table,
// optional groups, and zero padding. The output is vec(n x (m*q)).
real matrix rdistnn_native_loo(real matrix ranks, real matrix signal,
    real colvector table, | real colvector groups)
{
    string rowvector names
    real scalar n,m,q,grouped,length,count,rc,i
    real colvector payload
    real matrix result
    n=rows(ranks); m=cols(ranks); q=cols(signal); grouped=args()>=4
    if(n<2 | m<1 | q<1 | n>2147483647 | m>2147483647 | q>2147483647 |
        m*q>2147483647 | rows(signal)!=n | rows(table)!=n-1 |
        hasmissing(ranks) | hasmissing(signal) | hasmissing(table)) return(J(0,0,.))
    if(grouped) {
        if(rows(groups)!=n | hasmissing(groups)) return(J(0,0,.))
    }
    count=n*m*q
    length=5+n*m+n*q+n-1+grouped*n
    if(count>9007199254740991 | length>9007199254740991) return(J(0,0,.))
    payload=(1\n\m\q\grouped)\vec(ranks)\vec(signal)\table
    if(grouped) payload=payload\groups
    names=st_tempname(2)
    st_matrix(names[1],rdistnn_native_frame(payload))
    st_matrix(names[2],rdistnn_native_frame(J(count,1,.)))
    rc=_stata("plugin call rdistnn_plugin, loo_packed "+invtokens(names),1)
    result=J(0,0,.)
    if(rc==0) result=rdistnn_native_unframe(st_matrix(names[2]),n,m*q)
    for(i=1;i<=2;i++) st_matrix(names[i],J(0,0,.))
    if(rc) rdistnn_fail(rc,"Native deletion execution failed. Resolve the execution/resource error or select rdistnn_jackknife_engine mata explicitly.")
    if(rows(result)!=n | cols(result)!=m*q | hasmissing(result)) rdistnn_fail(498,"Native deletion execution returned malformed output. Reinstall the current plugin; if this persists, report a reproducible example.")
    return(result)
}
end
