# kt-3pi
Formalism to study J^PC -> 3 pi decays of mesons and rescattering effects.

Requires ROOT (tested with version 6.17) with the ROOT::Minuit2 library installed.

## Omega -> 3pi
The `param_fit` object allows an arbitrary amplitude to be fit to the a polynomial expansion around the center of the Dalitz plot. __scripts/omega_BW.cpp__ has an example with a Breit-Wigner amplitude used to extract plot parameters.   

Its the only script currently so to run simply do:
```
mkdir build && cd build
cmake ..
make && ./omega_BW
```
