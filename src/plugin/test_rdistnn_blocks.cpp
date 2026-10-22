#include "rdistnn_core.hpp"
#include <cassert>
#include <iostream>
#include <numeric>
#include <random>

static std::vector<double> calculate(const std::vector<int>& ranks,
    const std::vector<double>& signal, const std::vector<double>& table,
    const std::vector<int>& groups, int m=1, int q=1, bool direct=false)
{
    const int n=static_cast<int>(groups.size());
    std::vector<double> out(static_cast<std::size_t>(n)*m*q);
    if(direct) rdistnn_arithmetic::block_rank_loo(ranks.data(),signal.data(),table.data(),groups.data(),n,m,q,out.data(),true);
    else rdistnn::block_rank_loo(ranks.data(),signal.data(),table.data(),groups.data(),n,m,q,out.data());
    return out;
}

// Independent direct regrouping: sum represented weights in long double,
// round each actual averaged weight and product, then reduce products widely.
static std::vector<double> reference(const std::vector<int>& ranks,
    const std::vector<double>& signal, const std::vector<double>& table,
    const std::vector<int>& groups, int m, int q)
{
    const int n=static_cast<int>(groups.size());
    std::vector<double> out(static_cast<std::size_t>(n)*m*q);
    for(int target=0;target<m;++target) for(int deleted=0;deleted<n;++deleted) {
        std::vector<long double> mass(n,0);
        std::vector<int> count(n,0);
        for(int i=0;i<n;++i) if(i!=deleted) {
            int rank=ranks[i+target*n]-(ranks[i+target*n]>ranks[deleted+target*n]);
            mass[groups[i]-1]+=table[rank-1]; ++count[groups[i]-1];
        }
        for(int c=0;c<q;++c) {
            long double value=0;
            for(int i=0;i<n;++i) if(i!=deleted) {
                const int g=groups[i]-1;
                const double weight=static_cast<double>(mass[g]/count[g]);
                volatile double product=signal[i+c*n]*weight;
                value+=product;
            }
            out[deleted+(target*q+c)*n]=static_cast<double>(value);
        }
    }
    return out;
}

int main()
{
    const double large=std::ldexp(1.,54), eta=std::ldexp(1.,-1074);
    for(bool direct:{false,true}) {
        assert(calculate({1,2,3},{large,1,0},{1,0},{1,1,2},1,1,direct)==std::vector<double>({1,large,large/2}));
        assert(calculate({1,2,3,4},{large,1,large/2,large/2},{0,1,-1},{1,1,2,3},1,1,direct)==std::vector<double>({0,0,.5,.5}));
        assert(calculate({1,2,3},{eta,0,0},{1,0},{1,1,2},1,1,direct)==std::vector<double>({0,eta,0}));
        assert(calculate({1,2,3},{eta,0,0},{std::ldexp(1.,1000),0},{1,1,2},1,1,direct)==std::vector<double>({0,std::ldexp(1.,-74),std::ldexp(1.,-75)}));
        assert(calculate({1,2},{3,7},{-2},{1,2},1,1,direct)==std::vector<double>({-14,-6}));
    }
    const auto inter=calculate({1,2,3,4,5,6,7,8},{1,0,1,0,2,0,1,0},
        {large,0,1,0,0,-large,0},{1,2,1,3,1,4,1,5});
    assert(inter[4]==1);
    for(bool direct:{false,true}) {
        const auto observed=calculate({1,9,3,7,5,2,4,6,8,10},{0,1,1,1,1,0,0,0,0,0},
            {std::ldexp(1.,112),0,std::ldexp(1.,56),0,1,0,-std::ldexp(1.,56),0,-std::ldexp(1.,112)},
            {1,1,1,1,1,2,3,4,5,6},1,1,direct);
        assert(observed[9]==4./5.);
        const auto reordered=calculate({1,5,2,4,3,6},{0,1,1,1,1,0},
            {std::ldexp(1.,112),std::ldexp(1.,56),1,-std::ldexp(1.,56),-std::ldexp(1.,112)},
            {1,1,1,1,1,2},1,1,direct);
        assert(reordered[5]==4./5.);
    }
    // Every target/component layout, all-singletons, all-identical groups,
    // signed unnormalized tables and canonical sentinel boundaries.
    std::mt19937 generator(20261021);
    for(int trial=0;trial<200;++trial) {
        const int n=2+generator()%25,m=1+generator()%3,q=1+generator()%3;
        std::vector<int> groups(n),ranks(n*m),order(n);
        std::vector<double> signal(n*q),table(n-1);
        for(int i=0;i<n;++i) groups[i]=trial%3==0 ? 1 : (trial%3==1 ? i+1 : 1+i/3);
        for(double& v:signal) v=(static_cast<int>(generator()%201)-100)/4.;
        for(double& v:table) v=(static_cast<int>(generator()%21)-10)/16.;
        for(int target=0;target<m;++target) {
            std::iota(order.begin(),order.end(),0);
            if(trial%2) std::shuffle(order.begin(),order.end(),generator);
            else if(target%2) std::reverse(order.begin(),order.end());
            for(int r=0;r<n;++r) ranks[order[r]+target*n]=r+1;
        }
        const auto expected=reference(ranks,signal,table,groups,m,q);
        for(bool direct:{false,true}) {
            const auto actual=calculate(ranks,signal,table,groups,m,q,direct);
            for(std::size_t i=0;i<actual.size();++i)
                assert(std::abs(actual[i]-expected[i])<2e-12);
        }
    }
    bool threw=false;
    try { calculate({1,1,3},{1,2,3},{1,0},{1,1,2}); }
    catch(const std::invalid_argument&) { threw=true; }
    assert(threw);
    threw=false;
    try { calculate({1,2,3},{1,2,3},{1,0},{0,1,2}); }
    catch(const std::invalid_argument&) { threw=true; }
    assert(threw);
    for(int bad:{std::numeric_limits<int>::min(),std::numeric_limits<int>::max()}) {
        threw=false;
        try { calculate({bad,2,3},{1,2,3},{1,0},{1,1,2}); }
        catch(const std::invalid_argument&) { threw=true; }
        assert(threw);
    }
    std::cout << "Original-signal block deletion tests passed\n";
}
