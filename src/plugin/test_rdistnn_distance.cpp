#include "rdistnn_distance.hpp"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

namespace {
void require(bool ok, const char *what) {
    if (!ok) throw std::runtime_error(what);
}
template<class F> void invalid_unchanged(F fn, const std::vector<double>& out) {
    const auto before=out;
    bool rejected=false;
    try { fn(); } catch (const std::invalid_argument&) { rejected=true; }
    require(rejected,"invalid distance input accepted");
    require(out==before,"invalid distance input modified output");
}
// Independent wider-precision, target-first reference; no production helper use.
long double reference(const double* x,int n,int p,const double* target,int m,
                      int i,int q,const double* inv,bool l1) {
    long double total=0;
    for (int j=p-1;j>=0;--j) {
        long double d=static_cast<long double>(x[i+j*n])-target[q+j*m];
        if (inv) d*=inv[j];
        total+=l1?std::fabs(d):d*d;
    }
    return l1?total:std::sqrt(total);
}
}
int main(int argc, char** argv) {
    try {
        using rdistnn::distance_metric;
        const double x[]={0,3,-3, 0,4,4};
        const double targets[]={0,3,0, 0,0,0};
        const double scales[]={0.5,2.0};
        const double want_l1[]={0,7,7, 3,4,10, 0,7,7};
        const double want_l2[]={0,5,5, 3,4,std::sqrt(52.0), 0,5,5};
        std::vector<double> out(9),single(3);
        for (auto metric:{distance_metric::l1,distance_metric::l2}) {
            rdistnn::distance_batch_core(x,3,2,targets,3,nullptr,out.data(),metric);
            for (int j=0;j<9;++j)
                require(out[j]==(metric==distance_metric::l1?want_l1[j]:want_l2[j]),
                        "analytic unscaled/duplicate target distances");
            for (const double* inv:{static_cast<const double*>(nullptr),scales}) {
                rdistnn::distance_batch_core(x,3,2,targets,3,inv,out.data(),metric);
                for (int q=0;q<3;++q) {
                    const double target[]={targets[q],targets[q+3]};
                    rdistnn::distance_batch_core(x,3,2,target,1,inv,single.data(),metric);
                    for (int i=0;i<3;++i) {
                        const auto expected=reference(x,3,2,targets,3,i,q,inv,
                                                     metric==distance_metric::l1);
                        require(std::fabs(out[i+q*3]-expected)<1e-13,
                                "scaled analytic reference mismatch");
                        require(single[i]==out[i+q*3],"scalar/batch mismatch");
                    }

                }
            }
        }
        // Reordered targets must produce the corresponding column permutation.
        const double reordered[]={3,0,0, 0,0,0};
        rdistnn::distance_batch_core(x,3,2,reordered,3,nullptr,out.data(),distance_metric::l1);
        const std::vector<double> permuted={3,4,10,0,7,7,0,7,7};
        require(out==permuted,"target reorder changed distance values");
        // Aliasing is safe because all values are staged before publication.
        std::vector<double> aliased(x,x+6);
        const double alias_targets[]={0,3,0,0};
        rdistnn::distance_batch_core(aliased.data(),3,2,alias_targets,2,nullptr,
                                     aliased.data(),distance_metric::l1);
        require(aliased==std::vector<double>({0,7,7,3,4,10}),"aliased output changed inputs mid-computation");
        // Different matrix strides plus nontrivial scales catch transposition bugs.
        const int n=17,p=4,m=5;
        std::vector<double> rx(n*p),rt(m*p),rd(n*m);
        const double inverse[]={0.125,1.0,2.5,8.0};
        unsigned state=77123;
        for (auto* values:{&rx,&rt}) for (double& v:*values) {
            state=1664525u*state+1013904223u;
            v=(static_cast<int>((state>>8)%20001)-10000)/73.0;
        }
        for (auto metric:{distance_metric::l1,distance_metric::l2})
            for (const double* inv:{static_cast<const double*>(nullptr),inverse}) {
                rdistnn::distance_batch_core(rx.data(),n,p,rt.data(),m,inv,rd.data(),metric);
                for (int q=0;q<m;++q) for (int i=0;i<n;++i) {
                    const auto expected=reference(rx.data(),n,p,rt.data(),m,i,q,inv,
                                                 metric==distance_metric::l1);
                    require(std::fabs(rd[i+q*n]-expected)<=1e-12L*(1+expected),
                            "independent distance reference mismatch");
                }

            }
        // Distinct adjacent doubles remain distinct; exact distance ties stay equal.
        const double near[]={1,std::nextafter(1.0,2.0),-1,0};
        const double zero[]={0};
        double near_out[4];
        for (auto metric:{distance_metric::l1,distance_metric::l2}) {
            rdistnn::distance_batch_core(near,4,1,zero,1,nullptr,near_out,metric);
            require(near_out[0]==1 && near_out[2]==1 && near_out[3]==0 &&
                    near_out[1]==std::nextafter(1.0,2.0),"near tie distances changed");
        }
        out.assign(9,17.0);
        // FMA contraction can distinguish these coordinate-swapped rows even
        // though individually rounded products sum to an exact tie.
        // Build helpers pass these arguments so dimensions/data remain runtime
        // values: constant folding can conceal contraction in an FMA build.
        const int fma_n=argc==5?std::atoi(argv[1]):2;
        const int fma_p=argc==5?std::atoi(argv[2]):2;
        const double fma_a=argc==5?std::strtod(argv[3],nullptr):0.5790338377234229;
        const double fma_b=argc==5?std::strtod(argv[4],nullptr):1.4758492921057114;
        require(fma_n>=2 && fma_p>=2,"FMA fixture dimensions must be at least two");
        std::vector<double> fma_x(fma_n*fma_p),fma_zero(fma_p),fma_out(fma_n);
        fma_x[0]=fma_a; fma_x[1]=fma_b;
        fma_x[fma_n]=fma_b; fma_x[fma_n+1]=fma_a;
        rdistnn::distance_batch_core(fma_x.data(),fma_n,fma_p,fma_zero.data(),1,nullptr,
                                     fma_out.data(),distance_metric::l2);
        require(fma_out[0]==fma_out[1],"FMA changed a coordinate-swap distance tie");
        auto call=[&](const double* a,int nr,int nc,const double* b,int mr,
                      const double* inv,double* result,distance_metric metric=distance_metric::l2) {
            rdistnn::distance_batch_core(a,nr,nc,b,mr,inv,result,metric);
        };
        for (int bad:{0,-1}) {
            invalid_unchanged([&]{call(x,bad,2,targets,3,nullptr,out.data());},out);
            invalid_unchanged([&]{call(x,3,bad,targets,3,nullptr,out.data());},out);
            invalid_unchanged([&]{call(x,3,2,targets,bad,nullptr,out.data());},out);
        }
        invalid_unchanged([&]{call(nullptr,3,2,targets,3,nullptr,out.data());},out);
        invalid_unchanged([&]{call(x,3,2,nullptr,3,nullptr,out.data());},out);
        invalid_unchanged([&]{call(x,3,2,targets,3,nullptr,nullptr);},out);
        invalid_unchanged([&]{call(x,3,2,targets,3,nullptr,out.data(),static_cast<distance_metric>(99));},out);
        for (double bad:{0.0,-1.0,std::numeric_limits<double>::infinity(),
                         std::numeric_limits<double>::quiet_NaN()}) {
            const double bad_scales[]={1,bad};
            invalid_unchanged([&]{call(x,3,2,targets,3,bad_scales,out.data());},out);
        }
        for (double bad:{std::numeric_limits<double>::infinity(),
                         std::numeric_limits<double>::quiet_NaN()}) {
            double bad_x[]={0,3,-3,0,4,bad};
            invalid_unchanged([&]{call(bad_x,3,2,targets,3,nullptr,out.data());},out);
            invalid_unchanged([&]{call(x,3,2,bad_x,3,nullptr,out.data());},out);
        }
        invalid_unchanged([&]{call(x,std::numeric_limits<int>::max(),
                                   std::numeric_limits<int>::max(),targets,3,nullptr,out.data());},out);
        invalid_unchanged([&]{call(x,std::numeric_limits<int>::max(),1,
                                   targets,std::numeric_limits<int>::max(),nullptr,out.data());},out);
        const double extreme[]={0,1e308},negative[]={-1e308};
        invalid_unchanged([&]{call(extreme,2,1,negative,1,nullptr,out.data());},out);
        // Intermediate squares may overflow/underflow while the norm is finite.
        for(double h:{1e200,1e-200,std::ldexp(1.0,510),std::ldexp(1.0,-540)}) {
            const double huge[]={3*h,4*h,0,0,5*h, 4*h,3*h,h,0,0};
            const double z[]={0,0,0,0};
            const double expected[]={5,5,1,0,5};
            std::vector<double> distances(10), scalar(5);
            rdistnn::distance_batch_core(huge,5,2,z,2,nullptr,distances.data(),distance_metric::l2);
            rdistnn::distance_batch_core(huge,5,2,z,1,nullptr,scalar.data(),distance_metric::l2);
            for(int q=0;q<2;++q) for(int i=0;i<5;++i) {
                require(std::fabs(distances[i+q*5]/h-expected[i])<2e-14,"extreme analytic norm");
                require(distances[i+q*5]==scalar[i],"extreme scalar/batch parity");
            }
            double scaled[10];
            for(int i=0;i<10;++i) scaled[i]=huge[i]*(i<5?2:4);
            const double inverse[]={0.5,0.25};
            rdistnn::distance_batch_core(scaled,5,2,z,2,inverse,distances.data(),distance_metric::l2);
            for(int q=0;q<2;++q) for(int i=0;i<5;++i)
                require(distances[i+q*5]==scalar[i],"extreme inverse scaling");
        }
        // Stata reserves doubles at 2^1023 and above for missing values.
        // Every engine must choose the same exceptional arithmetic below it.
        const double boundary_x[]={7.385580696697091e153,1.076170170452885e154,
                                   7.827350851331151e153,0};
        const double boundary_zero[]={0,0};
        double boundary_out[2];
        call(boundary_x,2,2,boundary_zero,1,nullptr,boundary_out);
        require(boundary_out[0]==1.0761701704528852e154 &&
                boundary_out[1]==1.076170170452885e154,
                "common exceptional squared-sum boundary");
        const double boundary_scaled[]={8.000000000000001e307,8e307};
        const double boundary_target[]={-7.9e307},boundary_inverse[]={0.1};
        call(boundary_scaled,2,1,boundary_target,1,boundary_inverse,boundary_out);
        require(boundary_out[0]==1.5900000000000001e307 && boundary_out[1]==1.59e307,
                "common half-scaled subtraction preserves adjacent ranks");
        rdistnn::distance_batch_core(boundary_scaled,2,1,boundary_target,1,
            boundary_inverse,boundary_out,distance_metric::l1);
        require(boundary_out[0]==1.5900000000000001e307 && boundary_out[1]==1.59e307,
                "L1 common half-scaled subtraction preserves adjacent ranks");
        const double l1_large[]={1e308},l1_opposite[]={-1e308},l1_inverse[]={0.25};
        rdistnn::distance_batch_core(l1_large,1,1,l1_opposite,1,l1_inverse,
            boundary_out,distance_metric::l1);
        require(boundary_out[0]==5e307,"L1 scaled finite subtraction recovery");
        const double tiny=std::ldexp(1.0,-538);
        const std::vector<double> tiny_x(8,tiny),tiny_zero(8);
        rdistnn::distance_batch_core(tiny_x.data(),1,8,tiny_zero.data(),1,nullptr,out.data(),distance_metric::l2);
        require(out[0]==tiny*std::sqrt(8.0),"norm with subnormal squared sum lost");
        const double large[]={1e308},opposite[]={-1e308},small_inverse[]={0.25};
        call(large,1,1,opposite,1,small_inverse,out.data());
        require(out[0]==5e307,"scaled finite subtraction recovery");
        const double final_overflow[]={0,1.5e308,0,1.5e308},zero2[]={0,0};
        invalid_unchanged([&]{call(final_overflow,2,2,zero2,1,nullptr,out.data());},out);
        std::cout<<"rdistnn distance tests passed\n";
    } catch(const std::exception& e) {
        std::cerr<<"rdistnn distance test failed: "<<e.what()<<'\n'; return 1;
    }
}
