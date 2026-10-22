// Exercise the real packed transport entry point with independent SDK storage.
#include "rdistnn_plugin.cpp"
#include <iostream>
#include <map>

namespace {
struct matrix { int rows,cols; std::vector<double> values; };
std::map<std::string,matrix> matrices;
std::string failed_read;
int writes=0;
ST_int rowsof(char* name) {return matrices.count(name)?matrices.at(name).rows:0;}
ST_int colsof(char* name) {return matrices.count(name)?matrices.at(name).cols:0;}
ST_int readel(char* name,int r,int c,double* value) {
    if(failed_read==name) return 111;
    if(!matrices.count(name)) return 503;
    const auto& x=matrices.at(name);
    if(r<1 || c<1 || r>x.rows || c>x.cols) return 503;
    *value=x.values[(r-1)+static_cast<std::size_t>(c-1)*x.rows];return 0;
}
ST_int storeel(char* name,int r,int c,double value) {
    if(!matrices.count(name)) return 503;
    auto& x=matrices.at(name);
    if(r<1 || c<1 || r>x.rows || c>x.cols) return 503;
    ++writes;x.values[(r-1)+static_cast<std::size_t>(c-1)*x.rows]=value;return 0;
}
ST_int error(char*) {return 0;}
ST_boolean missing(double x) {return x>=8.9884656743115795e307;}
void require(bool value,const char* message) {if(!value) throw std::runtime_error(message);}
int invoke(std::initializer_list<const char*> args) {
    std::vector<std::string> strings;for(auto arg:args) strings.emplace_back(arg);
    std::vector<char*> argv;for(auto& arg:strings) argv.push_back(&arg[0]);
    return stata_call(static_cast<int>(argv.size()),argv.data());
}
matrix packed(const std::vector<double>& values) {
    const int rows=static_cast<int>(std::ceil(std::sqrt(values.size())));
    const int cols=static_cast<int>((values.size()+rows-1)/rows);
    matrix out={rows,cols,std::vector<double>(static_cast<std::size_t>(rows)*cols,0)};
    for(std::size_t k=0;k<values.size();++k)
        out.values[k/cols+(k%cols)*rows]=values[k];
    return out;
}
std::vector<double> unpack(const matrix& value,std::size_t length) {
    std::vector<double> out(length);
    for(std::size_t k=0;k<length;++k) out[k]=value.values[k/value.cols+(k%value.cols)*value.rows];
    return out;
}
void setcell(matrix& x,std::size_t k,double value) {x.values[k/x.cols+(k%x.cols)*x.rows]=value;}
void baseline(bool grouped=true) {
    matrices.clear();failed_read.clear();writes=0;
    const double large=std::ldexp(1.,54);
    std::vector<double> input={1,3,1,1,grouped?1.:0.,1,2,3,large,1,0,1,0};
    if(grouped) input.insert(input.end(),{1,1,2});
    matrices["F"]=packed(input);matrices["D"]=packed({17,17,17});
}
void rejected() {
    const auto before=matrices.at("D").values;const int count=writes;
    require(invoke({"loo_packed","F","D"})!=0,"invalid packed request accepted");
    require(writes==count && matrices.at("D").values==before,"invalid packed request modified output");
}
}
extern "C" {ST_plugin* _stata_=nullptr;}
int main() {
    ST_plugin api={};api.rowsof=rowsof;api.colsof=colsof;api.safematel=readel;
    api.safematstore=storeel;api.ismissing=missing;api.spouterr=error;_stata_=&api;
    try {
        baseline();matrices["API"]={1,1,{17}};
        require(invoke({"loo_packed_api","API"})==0,"packed native capability unavailable");
        require(matrices["API"].values[0]==1,"packed native capability version");
        require(invoke({"loo_packed_api","D"})!=0,"packed capability accepted invalid shape");
        require(invoke({"loo_packed","F","D"})==0,"packed grouped deletion failed");
        require(unpack(matrices["D"],3)==std::vector<double>({1,std::ldexp(1.,54),std::ldexp(1.,53)}),"packed grouped deletion values/layout");
        require(unpack(matrices["D"],4)[3]==0,"packed output padding not zeroed");
        // Existing raw operations and packed frames must agree bit for bit.
        for(bool grouped:{false,true}) {
            baseline(grouped);
            matrices["R"]={3,1,{1,2,3}};matrices["S"]={3,1,{std::ldexp(1.,54),1,0}};
            matrices["T"]={2,1,{1,0}};matrices["G"]={3,1,{1,1,2}};
            matrices["RAW"]={3,1,{17,17,17}};
            require((grouped?invoke({"block_rank_loo","R","S","T","G","RAW"}):invoke({"rank_loo","R","S","T","RAW"}))==0,"raw comparison failed");
            require(invoke({"loo_packed","F","D"})==0,"packed comparison failed");
            require(unpack(matrices["D"],3)==matrices["RAW"].values,"raw/packed arithmetic differs");
        }
        // Multiple targets and components retain the raw column layout.
        for(bool grouped:{false,true}) for(int n:{6,17}) {
            baseline(grouped);
            const int m=2,q=3;
            std::vector<double> ranks,signal,table,groups;
            for(int j=0;j<m;++j) for(int i=0;i<n;++i) ranks.push_back(j?n-i:i+1);
            for(int j=0;j<q;++j) for(int i=0;i<n;++i) signal.push_back(std::sin(i+j*3.0));
            for(int i=0;i<n-1;++i) table.push_back((i%3-1)*0.25);
            for(int i=0;i<n;++i) groups.push_back(i/3+1);
            std::vector<double> input={1.,static_cast<double>(n),2.,3.,grouped?1.:0.};
            input.insert(input.end(),ranks.begin(),ranks.end());
            input.insert(input.end(),signal.begin(),signal.end());
            input.insert(input.end(),table.begin(),table.end());
            if(grouped) input.insert(input.end(),groups.begin(),groups.end());
            matrices["F"]=packed(input);matrices["D"]=packed(std::vector<double>(n*m*q,17));
            matrices["R"]={n,m,ranks};matrices["S"]={n,q,signal};matrices["T"]={n-1,1,table};
            matrices["G"]={n,1,groups};matrices["RAW"]={n,m*q,std::vector<double>(n*m*q,17)};
            require((grouped?invoke({"block_rank_loo","R","S","T","G","RAW"}):invoke({"rank_loo","R","S","T","RAW"}))==0,"raw multi-dimensional comparison failed");
            require(invoke({"loo_packed","F","D"})==0,"packed multi-dimensional comparison failed");
            require(unpack(matrices["D"],n*m*q)==matrices["RAW"].values,"raw/packed multi-dimensional layout differs");
        }
        baseline();
        matrices["F"]=packed({1,4,1,1,1,1,2,3,4,1,2,3,4,1,0,0,1,2,1,2});
        matrices["D"]=packed({17,17,17,17});
        require(invoke({"loo_packed","F","D"})==0,"packed interleaved groups failed");
        require(unpack(matrices["D"],4)==std::vector<double>({3,2,1,2}),"packed interleaved survivor values");
        // An input/output alias with identical balanced shapes is safe.
        baseline(false);
        matrices["F"]=packed({1,2,2,5,0,1,2,2,1,1,2,3,4,5,6,7,8,9,10,1});
        require(invoke({"loo_packed","F","F"})==0,"packed input/output alias failed");
        require(unpack(matrices["F"],20)==std::vector<double>({2,1,4,3,6,5,8,7,10,9,2,1,4,3,6,5,8,7,10,9}),"packed alias values/layout");
        for(std::size_t field=0;field<5;++field) for(double bad:{-1.,.5,std::numeric_limits<double>::quiet_NaN(),std::numeric_limits<double>::infinity()}) {
            baseline();setcell(matrices["F"],field,bad);rejected();
        }
        for(double bad:{0.,1.,static_cast<double>(INT_MAX)+1.,std::numeric_limits<double>::max()}) {
            baseline();setcell(matrices["F"],1,bad);rejected();
        }
        baseline();setcell(matrices["F"],2,INT_MAX);setcell(matrices["F"],3,INT_MAX);rejected();
        baseline();setcell(matrices["F"],1,INT_MAX);setcell(matrices["F"],2,10000000);rejected();
        for(double bad:{0.,-2147483648.,1.5,4.,std::numeric_limits<double>::quiet_NaN()}) {
            baseline();setcell(matrices["F"],5,bad);rejected();
        }
        baseline();setcell(matrices["F"],6,1);rejected();
        baseline();setcell(matrices["F"],8,std::numeric_limits<double>::quiet_NaN());rejected();
        baseline();setcell(matrices["F"],13,0);rejected();
        baseline(false);setcell(matrices["F"],15,1);rejected();
        baseline();matrices["F"].rows=2;matrices["F"].cols=8;rejected();
        baseline();matrices["D"]={1,3,{17,17,17}};rejected();
        baseline();failed_read="F";const int count=writes;
        require(invoke({"loo_packed","F","D"})==111 && writes==count,"packed read failure not propagated");
        std::cout<<"Packed native deletion protocol tests passed\n";
    } catch(const std::exception& error) {std::cerr<<error.what()<<'\n';return 1;}
}
