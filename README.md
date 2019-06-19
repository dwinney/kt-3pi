# kt_3pi
Codebase with tools to study J^PC -> 3 pi decays of mesons and rescattering effects with Khuri-Treiman formalism.

Requires ROOT (tested with version 6.17) with the ROOT::Minuit2 and MathMore libraries installed.
MathMore additionally requires a distribution of GSL (version >= 1.8) linked with ROOT.

## Dalitz Plot Parameter Extraction
The `param_fit` object allows an arbitrary amplitude to be fit to the a polynomial expansion around the center of the Dalitz plot. __omega_BW.cpp__ has an example with a Breit-Wigner amplitude used to extract plot parameters for the omega meson decaying to three pions.   

To build in command-line use:
```
mkdir build && cd build
cmake ..
make omega_BW && ./omega_BW
```

Additionally you can make your own amplitude and use the `param_fit` and `dalitz` objects to extract parameters. Any object can be used so long as has a `decay_kinematics` object named "kinematics" as a publicly accessible data member and is callable with the right signature (```double ()(double, double)```
or
```std::complex<double> ()(double, double)```)
to evaluate the user-defined amplitude at some values of Mandelstam s and t. Then an 'output_amp' object with 'num_param' free fitting parameters is output with:

```
my_amp example(my_inputs);
dalitz_fit<my_amp, output_amp> fit(example);

output_amp = fit.extract_params(num_params);
```

## Khuri-Treiman Equations
__UNDER CONSTRUCTION__

The ```kt_amplitude``` class will allow to create dispersive isobar amplitudes within the KT formalism for the three pion final state.
The class structure will be laid out as:

Mass and quantum numbers of decaying particle which define the amplitude ( ``` kt_amplitude``` )

 -> Number of partial waves to consider (each with a spin and isospin projection to be summed over) ( ``` isobar``` )

 -> Number of subtractions in the dispersion relation and iterations of the rescattering ladder ( ```iteration``` )

 -> Method of evaluating the KT-equations ( currently only ```conformal_int``` which uses the method in [[Dan14a]](https://arxiv.org/abs/1409.7708) but a Pasquier-inversion method of [[Guo14]](https://arxiv.org/abs/1412.3970) can easily be implemented).

 #### Example
The __omega_KT.cpp__ script gives an example calculation of the _unsubtracted_ KT equations for the omega meson case for arbitrary iterations of the integral equation. To build use see the example above with __omega_KT__ instead of __omega_BW__.

Because there is only one relevent isobar and no isospin dependence, we use the `isobar` class directly but for more complicated cases we will use a combination of `isobar` objects in the `kt_amplitude` class.

Subtracted dispersion relations will soon be fittable to data using the above `dalitz_fit` class.
