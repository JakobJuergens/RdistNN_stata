// Independent analytic checks for asymptotic Beta interval rank weights.
#include "rdistnn_core.hpp"
#include <iostream>
#include <numeric>

namespace {
int failures=0;
void check(bool value,const char* label) {
    if (!value) { std::cerr << "FAIL: " << label << '\n'; ++failures; }
}
bool near_relative(double got,long double expected,double tolerance) {
    return std::isfinite(got) && got>0 &&
        std::fabs(static_cast<long double>(got)/expected-1)<tolerance;
}
double factorial(int n) {
    double value=1;
    for(int i=2;i<=n;++i) value*=i;
    return value;
}
double choose_integer(int n,int k) {
    return factorial(n)/(factorial(k)*factorial(n-k));
}
// Independently integrate the integer-shape Beta density on a local interval.
// Set x=lower+h*t and 1-x=(1-upper)+h*(1-t), expand both integer powers,
// then integrate t^i*(1-t)^j as i!*j!/(i+j+1)!. Every term is nonnegative.
// Deliberately use double: MSVC has no wider long double, and subtracting
// two near-one binomial CDFs would invalidate a small-tail reference there.
double beta_interval_polynomial(int a,int b,double lower,double upper) {
    const double h=upper-lower;
    double integral=0;
    for(int i=0;i<a;++i) {
        for(int j=0;j<b;++j) {
            integral+=choose_integer(a-1,i)*choose_integer(b-1,j)
                *std::pow(lower,a-1-i)*std::pow(1-upper,b-1-j)
                *std::pow(h,i+j)*factorial(i)*factorial(j)/factorial(i+j+1);
        }
    }
    return integral*h*factorial(a+b-1)/(factorial(a-1)*factorial(b-1));
}
}
int main() {
    for(int n=2;n<=8;++n) for(int s=1;s<=n;++s) for(int k=1;k<=s;++k) {
        const auto table=rdistnn_weights::exact_table(n,k,s);
        for(int r=1;r<=n;++r) {
            double expected=0;
            for(int l=1;l<=k;++l) if(l<=r && s-l<=n-r)
                expected+=choose_integer(r-1,l-1)*choose_integer(n-r,s-l)/choose_integer(n,s)/k;
            check(std::fabs(table[r-1]-expected)<3e-15,"exact order statistic combinatorial oracle");
        }
    }
    const auto exact_tail=rdistnn_weights::exact_table(2000,1,500);
    for(int i=1;i<2000;++i) check(exact_tail[i]<=exact_tail[i-1],"scaled exact tail monotonicity");

    // The fixed-budget continued fraction cannot resolve this central tail.
    // Returning its last iterate would silently violate the exact symmetry
    // I_0.5(a,a)=0.5. Exhaustion must propagate through interval evaluation.
    const auto beta_failure=[](auto evaluate) {
        try { evaluate(); }
        catch(const std::runtime_error& error) {
            return std::string(error.what()).find("Beta")!=std::string::npos;
        }
        return false;
    };
    for(bool upper:{false,true}) {
        check(beta_failure([upper] { rdistnn::beta_probability(1e7,1e7,.5,upper); }),
              "nonconverged Beta tail reports a numerical failure");
    }
    check(beta_failure([] { rdistnn::beta_interval_probability(1e7,1e7,.49,.5); }),
          "Beta interval propagates a failed endpoint evaluation");
    check(std::fabs(rdistnn::beta_probability(1000,1000,.5)-.5)<1e-10,
          "converged symmetric Beta remains available");
    check(rdistnn::beta_probability(1e7,1e7,0)==0 &&
          rdistnn::beta_probability(1e7,1e7,1)==1 &&
          rdistnn::beta_probability(1e7,1e7,0,true)==1 &&
          rdistnn::beta_probability(1e7,1e7,1,true)==0,
          "large-shape exact endpoints bypass iteration");
    check(rdistnn::beta_probability(1,1000,.99,true)==0,
          "legitimate tail underflow remains zero");
    const int n=100,s=80;
    std::vector<int> ranks(n); std::iota(ranks.begin(),ranks.end(),1);
    std::vector<double> weights(n);
    rdistnn::dnn_weights_asymp_from_ranks(ranks.data(),n,1,s,weights.data());
    for(int r : {50,80,99,100}) {
        // Use the exact double endpoints consumed by the production kernel.
        const long double lo=static_cast<double>(r-1)/n,hi=static_cast<double>(r)/n;
        const long double expected=std::pow(1-lo,s)-std::pow(1-hi,s);
        check(near_relative(weights[r-1],expected,2e-12),"Beta(1,80) tiny tail survives");
    }
    check(std::fabs(std::accumulate(weights.begin(),weights.end(),0.0)-1)<2e-14,
          "tail-preserving rank table normalizes");
    std::vector<double> X(n); for(int i=0;i<n;++i) X[i]=i;
    double x0=0;
    rdistnn::dnn_weights_l2(X.data(),n,1,&x0,nullptr,1,s,weights.data(),true);
    check(near_relative(weights[99],std::pow(1-static_cast<long double>(99.0/100),80),2e-12),
          "legacy weights_l2 retains final tiny mass");
    // Both small tails, crossing the mean, zero-width bins, and endpoints.
    for(double lo : {0.0,.2,.79,.99}) {
        const double hi=lo==.99 ? 1.0 : lo+.01;
        const long double expected=std::pow(1-static_cast<long double>(lo),100)
            -std::pow(1-static_cast<long double>(hi),100);
        check(near_relative(rdistnn::beta_interval_probability(1,100,lo,hi),expected,3e-12),
              "direct Beta(1,100) analytic tail");
        // Reflect exactly the double endpoints passed to the helper.
        const double rlo=1-hi,rhi=1-lo;
        const long double reflected_expected=std::pow(static_cast<long double>(rhi),100)
            -std::pow(static_cast<long double>(rlo),100);
        check(near_relative(rdistnn::beta_interval_probability(100,1,rlo,rhi),reflected_expected,3e-12),
              "reflected Beta(100,1) analytic tail");
    }
    for(int a=1;a<=7;++a) {
        check(rdistnn::beta_interval_probability(a,8-a,0,1)==1,
              "full unit interval has exactly unit mass");
        check(rdistnn::beta_interval_probability(a,8-a,0,0)==0 &&
              rdistnn::beta_interval_probability(a,8-a,1,1)==0,
              "zero-width endpoint intervals have zero mass");
        double total=0;
        for(int r=1;r<=24;++r) total+=rdistnn::beta_interval_probability(a,8-a,(r-1)/24.0,r/24.0);
        check(std::fabs(total-1)<2e-14,"every integer shape table sums to one");
    }
    const double subnormal=rdistnn::beta_interval_probability(1,160,.99,1);
    const double subnormal_ref=static_cast<double>(std::pow(1-static_cast<long double>(.99),160));
    check(subnormal>0 && std::fabs(subnormal-subnormal_ref)<=8*std::numeric_limits<double>::denorm_min(),
          "representable subnormal tail is never cut off");
    const int ordinary_n=24,ordinary_s=7;
    ranks.resize(ordinary_n); std::iota(ranks.rbegin(),ranks.rend(),1);
    weights.resize(ordinary_n);
    for(int k=1;k<=ordinary_s;++k) {
        rdistnn::dnn_weights_asymp_from_ranks(ranks.data(),ordinary_n,k,ordinary_s,weights.data());
        for(int i=0;i<ordinary_n;++i) {
            int r=ranks[i];
            const long double lo=static_cast<double>(r-1)/ordinary_n;
            const long double hi=static_cast<double>(r)/ordinary_n;
            long double expected=0;
            for(int a=1;a<=k;++a) expected+=beta_interval_polynomial(a,ordinary_s-a+1,static_cast<double>(lo),static_cast<double>(hi))/k;
            check(near_relative(weights[i],expected,5e-11),"all k, reversed ranks, independent polynomial integral oracle");
        }
    }
    if(failures) return 1;
    std::cout << "Stable Beta weight tests passed.\n";
}
