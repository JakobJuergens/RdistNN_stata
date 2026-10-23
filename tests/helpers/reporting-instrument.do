version 17
// Test-only source instrumentation: the real public callers and bridge execute.
// No production counter/global or native testing operation is installed.
capture mata: mata drop rdistnn_test_instrument()
mata:
void rdistnn_test_instrument(string scalar source, string scalar destination,
    string scalar before, string scalar after)
{
    real scalar input,output
    string scalar line
    input=fopen(source,"r"); unlink(destination); output=fopen(destination,"w")
    while((line=fget(input))!=J(0,0,"")) fput(output,subinstr(line,before,after))
    fclose(input); fclose(output)
}
end

capture program drop _rdistnn_test_notice_count
program _rdistnn_test_notice_count, rclass
    syntax using/, PATTERN(string)
    tempname handle
    file open `handle' using "`using'", read text
    local count 0
    file read `handle' line
    while r(eof)==0 {
        if strpos(`"`line'"',`"`pattern'"') local ++count
        file read `handle' line
    }
    file close `handle'
    return scalar count=`count'
end
