#include "omnes.hpp"

int main()
{
        omnes pwave(1,1);
        pwave.set_N_omnes(50);
        pwave.set_eps(0, 0);
        std::cout << pwave.eval(1.) << "\n";
        return 0;
}
