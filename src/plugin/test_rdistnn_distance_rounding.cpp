// Keep this small translation unit separate: runtime dimensions and data
// expose contraction that constant-folded analytic fixtures can conceal.
#include "rdistnn_distance.hpp"
#include <cstdlib>
#include <iostream>
#include <vector>

int main(int argc, char** argv)
{
    if(argc!=5) return 2;
    const int n=std::atoi(argv[1]), p=std::atoi(argv[2]);
    if(n<2 || p<2 || n>100 || p>100) return 2;
    const double a=std::strtod(argv[3],nullptr), b=std::strtod(argv[4],nullptr);
    std::vector<double> x(n*p),target(p),out(n);
    x[0]=a; x[1]=b; x[n]=b; x[n+1]=a;
    rdistnn::distance_batch_core(x.data(),n,p,target.data(),1,nullptr,
                                 out.data(),rdistnn::distance_metric::l2);
    if(out[0]!=out[1]) {
        std::cerr << "L2 contraction changed a coordinate-swap tie\n";
        return 1;
    }
    std::cout << "rdistnn distance rounding test passed\n";
}
