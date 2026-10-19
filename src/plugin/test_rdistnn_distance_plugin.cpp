// Execute the real entry point against a narrow host SDK callback fixture.
// Real Stata integration remains covered by the package's runtime .do tests.
#include "rdistnn_core.hpp"
namespace { bool probe_beta_failure=false; }
namespace rdistnn {
// Keep the real Beta failure and real plugin boundary. A full rank-weight
// table reaching these shapes would require millions of rows and a huge k.
// Substitute that expensive call only in this host-protocol test.
void weights_l2_beta_probe(const double* X,int n,int p,const double* x0,
    const double* inv_s,int k,int s,double* out,bool asymp) {
    if(probe_beta_failure) beta_probability(1e7,1e7,.5);
    dnn_weights_l2(X,n,p,x0,inv_s,k,s,out,asymp);
}
}
#define dnn_weights_l2 weights_l2_beta_probe
#include "rdistnn_plugin.cpp"
#undef dnn_weights_l2
#include <iostream>
#include <map>
#include <stdexcept>
#include <utility>

namespace {
struct matrix { int n,p; std::vector<double> values; };
std::map<std::string,matrix> matrices;
std::string read_failure;
std::string last_error;
std::vector<double> data_x,data_w;
int data_writes=0;
int reads=0,writes=0;
ST_int nobs() { return 3; }
ST_int nvars() { return 2; }
ST_int in1() { return 1; }
ST_int in2() { return 3; }
ST_boolean selected(ST_int) { return 1; }
ST_int readvar(ST_int j,ST_int i,double* value) {
    if(i<1 || i>3 || j<1 || j>2) return 503;
    *value=(j==1 ? data_x : data_w).at(i-1); return 0;
}
ST_int storevar(ST_int j,ST_int i,double value) {
    if(j!=2 || i<1 || i>3) return 503;
    ++data_writes; data_w.at(i-1)=value; return 0;
}
ST_int rowsof(char* name) { auto m=matrices.find(name); return m==matrices.end()?0:m->second.n; }
ST_int colsof(char* name) { auto m=matrices.find(name); return m==matrices.end()?0:m->second.p; }
ST_int readel(char* name,int i,int j,double* v) {
    ++reads;
    if (read_failure==name) return 111;
    const auto m=matrices.find(name);
    if (m==matrices.end() || i<1 || j<1 || i>m->second.n || j>m->second.p) return 503;
    *v=m->second.values.at((i-1)+static_cast<std::size_t>(j-1)*m->second.n); return 0;
}
ST_int storeel(char* name,int i,int j,double v) {
    ++writes;
    auto m=matrices.find(name);
    if (m==matrices.end() || i<1 || j<1 || i>m->second.n || j>m->second.p) return 503;
    m->second.values.at((i-1)+static_cast<std::size_t>(j-1)*m->second.n)=v; return 0;
}
ST_int error(char* message) {last_error=message; return 0;}
ST_boolean missing(double v) {return v>=8e307;}
void require(bool ok,const char* msg) { if(!ok) throw std::runtime_error(msg); }
int invoke(std::initializer_list<const char*> args) {
    std::vector<std::string> strings; for(auto a:args) strings.emplace_back(a);
    std::vector<char*> argv; for(auto& s:strings) argv.push_back(&s[0]);
    return stata_call(static_cast<int>(argv.size()),argv.data());
}
void reset() {
    matrices.clear(); read_failure.clear(); reads=0; writes=0;
    last_error.clear(); probe_beta_failure=false;
    data_x={0,1,2}; data_w={73,73,73}; data_writes=0;
    matrices["X"]={3,2,{0,3,-3,0,4,4}};
    matrices["X0"]={2,2,{0,3,0,0}};
    matrices["S"]={1,2,{0.5,2}};
    matrices["D"]={3,2,std::vector<double>(6,17)};
}
void rejected(std::initializer_list<const char*> args) {
    const auto before=matrices.at("D").values; const int before_writes=writes;
    require(invoke(args)!=0,"malformed plugin request accepted");
    require(matrices.at("D").values==before && writes==before_writes,
            "rejected plugin request published output");
}
}
extern "C" { ST_plugin* _stata_=nullptr; }
int main() {
    ST_plugin api={};
    api.rowsof=rowsof; api.colsof=colsof; api.safematel=readel;
    api.safematstore=storeel; api.ismissing=missing; api.spouterr=error;
    api.nobs=nobs; api.nvars=nvars; api.nobs1=in1; api.nobs2=in2;
    api.selobs=selected; api.safevdata=readvar; api.safestore=storevar;
    _stata_=&api;
    try {
        reset(); matrices["V"]={1,1,{17}};
        require(invoke({"distance_api","V"})==0,"distance capability operation unavailable");
        require(matrices["V"].values[0]==1,"wrong distance capability version");
        require(invoke({"distance_api","D"})!=0,"capability probe accepted wrong output shape");
        require(matrices["D"].values==std::vector<double>(6,17),"capability probe modified invalid output");
        reset();
        require(invoke({"distance_batch","l1","X","X0","none","D"})==0,"L1 protocol failed");
        require(matrices["D"].values==std::vector<double>({0,7,7,3,4,10}),"L1 protocol values/layout");
        require(invoke({"distance_batch","l2","X","X0","S","D"})==0,"scaled L2 protocol failed");
        const std::vector<double> want={0,std::sqrt(66.25),std::sqrt(66.25),1.5,8,std::sqrt(73.0)};
        require(matrices["D"].values==want,"scaled L2 protocol values/layout");
        reset(); matrices["X0"]={1,2,{0,0}}; matrices["D"]={3,1,{17,17,17}};
        require(invoke({"distance_batch","l2","X","X0","none","D"})==0,"scalar target protocol failed");
        require(matrices["D"].values==std::vector<double>({0,5,5}),"scalar target protocol values");
        reset(); read_failure="S";
        require(invoke({"distance_batch","l1","X","X0","S","D"})==111 && writes==0,
                "inverse-scale read failure published output or lost return code");
        reset();
        rejected({"distance_batch","l1","X","X0","D"});
        rejected({"distance_batch","bad","X","X0","none","D"});
        rejected({"distance_batch","l2","absent","X0","none","D"});
        matrices["X0"].p=1;
        rejected({"distance_batch","l2","X","X0","none","D"}); reset();
        matrices["S"]={2,1,{1,1}};
        rejected({"distance_batch","l2","X","X0","S","D"}); reset();
        for(double bad:{0.0,-1.0,8.988e307,std::numeric_limits<double>::infinity(),
                        std::numeric_limits<double>::quiet_NaN()}) {
            reset(); matrices["S"].values[1]=bad;
            rejected({"distance_batch","l2","X","X0","S","D"});
        }
        for(const char* which:{"X","X0"}) {
            reset(); matrices[which].values.back()=8.988e307;
            rejected({"distance_batch","l2","X","X0","none","D"});
            reset(); read_failure=which;
            require(invoke({"distance_batch","l2","X","X0","none","D"})==111,
                    "matrix read error not propagated");
            require(writes==0,"read failure published output");
        }
        reset(); matrices["D"]={2,3,std::vector<double>(6,17)};
        rejected({"distance_batch","l1","X","X0","none","D"});
        require(reads==0,"wrong output shape read inputs before validation");
        reset(); matrices["X"]={2,1,{0,7e307}}; matrices["X0"]={1,1,{-7e307}};
        matrices["D"]={2,1,{17,17}};
        rejected({"distance_batch","l1","X","X0","none","D"});
        for(double h:{1e200,1e-200}) {
            reset(); matrices["X"]={2,2,{3*h,4*h,4*h,3*h}};
            matrices["X0"]={1,2,{0,0}}; matrices["D"]={2,1,{17,17}};
            require(invoke({"distance_batch","l2","X","X0","none","D"})==0,
                    "representable extreme L2 protocol failed");
            for(double d:matrices["D"].values)
                require(std::fabs(d/h-5)<2e-14,"extreme L2 analytic protocol values");
        }
        reset(); matrices["X"]={2,2,{0,7e307,0,7e307}};
        matrices["X0"]={1,2,{0,0}}; matrices["D"]={2,1,{17,17}};
        rejected({"distance_batch","l2","X","X0","none","D"});
        reset(); matrices["X"].n=std::numeric_limits<int>::max();
        matrices["X"].p=std::numeric_limits<int>::max();
        matrices["X0"].p=std::numeric_limits<int>::max();
        matrices["D"].n=std::numeric_limits<int>::max();
        rejected({"distance_batch","l2","X","X0","none","D"});
        require(reads==0,"oversized product read inputs");
        require(stata_call(0,nullptr)!=0,"empty entry point accepted");
        reset();
        require(invoke({"weights_l2","asymp","1","2","0"})==0,
                "converged native Beta weights failed");
        const double beta_weights[3]={5./9,1./3,1./9};
        for(int i=0;i<3;++i)
            require(std::fabs(data_w[i]-beta_weights[i])<1e-14,
                    "converged native Beta weights changed");
        reset(); probe_beta_failure=true;
        require(invoke({"weights_l2","asymp","1","2","0"})==498,
                "Beta nonconvergence did not report Stata numerical error 498");
        require(data_writes==0 && data_w==std::vector<double>({73,73,73}),
                "failed Beta calculation modified caller output");
        require(last_error.find("Beta")!=std::string::npos &&
                last_error.find("Mata")!=std::string::npos,
                "Beta failure omitted its context or recovery guidance");
        std::cout<<"rdistnn distance plugin protocol tests passed\n";
    } catch(const std::exception& e) {
        std::cerr<<"rdistnn plugin test failed: "<<e.what()<<'\n'; return 1;
    }
}
