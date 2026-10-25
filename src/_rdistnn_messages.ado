*! version 0.1.5  2026-10-25
version 17
// Native Stata channels. No r()/e() diagnostic store and no session-wide scope.
capture program drop _rdistnn_messages
program define _rdistnn_messages
    version 17
    syntax anything(name=kind), MESSAGE(string) [HINT(string) RC(integer 498)]
    if !inlist("`kind'","error","warning","note") exit 198
    if "`kind'"!="error" & !c(noisily) exit
    // Normalize ASCII controls without loading Mata or changing r()/e().
    local message = ustrregexra(`"`macval(message)'"', "[\x00-\x1f]", " ")
    local hint = ustrregexra(`"`macval(hint)'"', "[\x00-\x1f]", " ")
    // _asis keeps dynamic braces, percent signs and quotes literal.
    if "`kind'"=="error" {
        display as error _asis `"`message' `hint'"'
        exit `rc'
    }
    display as text _asis `"`kind': `message' `hint'"'
end

capture mata: mata drop rdistnn_fail()
capture mata: mata drop rdistnn_message_value()
mata:
string scalar rdistnn_message_value(string scalar value)
{
    real scalar i
    for(i=0;i<32;i++) value=subinstr(value,char(i)," ")
    if(ustrlen(value)>80) value=usubstr(value,1,77)+"..."
    return(value)
}

void rdistnn_fail(real scalar rc, string scalar message)
{
    real scalar i
    for(i=0;i<32;i++) message=subinstr(message,char(i)," ")
    displayas("error")
    display(message,1)
    exit(rc)
}
end
