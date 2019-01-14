#include "omnes.hpp"
#include <iomanip>

int main()
{
        omnes pwave(1,1);
        pwave.set_N_omnes(400);

        for (int i; i < 100; i++)
        {
                std::cout << std::left << std::setw(10) << sqrt(0.01 *double(i)) << std::setw(40) << pwave.eval(0.01*double(i)) << "\n";
        }

        return 0;
}
