#include "omnes.hpp"
#include <iomanip>

// std::complex<double> test(double x, int e)
// {
//   return x ;
// }

int main()
{
        // amplitude a;
        // pipi one(1);
        // std::cout << amplitude::sthPi<< "\n";

        // omnes oneone(1,1);
        omnes pwave(1,1);
        pwave.set_N_omnes(400);
        
        for (int i; i < 100; i++)
        {
                std::cout << std::left << std::setw(10) << sqrt(0.01 *double(i)) << std::setw(40) << pwave.eval(0.01*double(i)) << "\n";
        }

        // std::complex<double> (*pt_test)(double, int);
        // pt_test = test;
        //
        // kt_amplitude my_kt(test);
        // std::cout << my_kt.eval(5., 2) << "\n";
        // my_kt.start();

        return 0;
}
