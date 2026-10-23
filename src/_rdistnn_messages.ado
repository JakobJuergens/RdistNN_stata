version 17
// Native Stata channels. No r()/e() diagnostic store and no session-wide scope.
capture program drop _rdistnn_messages
program define _rdistnn_messages
    version 17
    syntax anything(name=kind), MESSAGE(string) [HINT(string) RC(integer 498)]
    if !inlist("`kind'","error","warning","note") exit 198
    if "`kind'"!="error" & !c(noisily) exit
    // _asis keeps dynamic braces, percent signs and quotes literal.
    if "`kind'"=="error" {
        display as error _asis `"`message' `hint'"'
        exit `rc'
    }
    display as text _asis `"`kind': `message' `hint'"'
end

capture mata: mata drop rdistnn_fail()
capture mata: mata drop rdistnn_message_value()
capture mata: mata drop rdistnn_message_examples()
mata:
string scalar rdistnn_message_value(string scalar value)
{
    real scalar i
    for(i=0;i<32;i++) value=subinstr(value,char(i)," ")
    if(ustrlen(value)>80) value=usubstr(value,1,77)+"..."
    return(value)
}
string scalar rdistnn_message_examples(real colvector ids)
{
    string scalar out
    real scalar i,n
    n=rows(ids); out=""
    for(i=1;i<=min((5,n));i++) out=out+(i>1 ? ", " : "")+strofreal(ids[i],"%12.0g")
    if(n>5) out=out+" ("+strofreal(n-5,"%12.0g")+" more)"
    return(out)
}
void rdistnn_fail(real scalar rc, string scalar message)
{
    displayas("error")
    display(message)
    exit(rc)
}
end
