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

Additionally you can make your own amplitude and use the `param_fit` and `dalitz` objects to extract parameters. Any object can be used so long as it inherits publicly from the `amplitude` class and is callable with the right signature (```double ()(double, double)```
or
```std::complex<double> ()(double, double)```)
to evaluate the user-defined amplitude at some values of Mandelstam s and t. Then an 'output_amp' object with 'num_param' free fitting parameters is output with:

```
my_amp example(my_inputs);
dalitz_fit<my_amp, output_amp> fit(example);

output_amp = fit.extract_params(num_params);
```
