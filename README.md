# kt_3pi
Codebase with tools to study J^PC -> 3 pi decays of mesons and rescattering effects with Khuri-Treiman formalism.

Requires ROOT (tested with version 6.17) with the ROOT::Minuit2 and MathMore libraries installed.
MathMore additionally requires a distribution of GSL (version >= 1.8) linked with ROOT.

## Dalitz Plot Parameter Extraction
The `dalitz_fit` object allows an arbitrary amplitude to be fit to the a polynomial expansion around the center of the Dalitz plot. __omega_BW.cpp__ contains an example with a Breit-Wigner amplitude used to extract plot parameters for the omega meson decaying to three pions.   
To build in command-line use:
```
mkdir build && cd build
cmake ..
make omega_BW && ./omega_BW
```

Additionally you can make your own amplitude and use the `dalitz` and `dalitz_fit` objects to extract parameters. Any object can be used so long as has the following publicly-accessible members:
```
decay_kinematics kinematics;
void set_params(n_params, const * par);
complex<double> eval(double s, double t);
```
To store kinematic quantities of the reaction, appropriately set the fitting parameters from an input array, and evaluate the amplitude at Mandelstam s and t in the decay region  respectively .
```
my_amp example(my_inputs);
output_amp results(my_inputs);

dalitz_fit<my_amp, output_amp> fit(&example, &results);
fit.extract_params(num_params);
```

## Khuri-Treiman Equations
__UNDER CONSTRUCTION__

The ```kt_amplitude``` class will allow to create dispersive isobar amplitudes within the KT formalism for the three pion final state.
The class structure will be laid out as:

Mass and quantum numbers of decaying particle which define the amplitude ( ``` kt_amplitude``` )

 -> Number of partial waves to consider (each with a spin and isospin projection to be summed over) ( ``` isobar``` )

 -> Number of subtractions in the dispersion relation and iterations of the rescattering ladder ( ```iteration``` )

 -> Method of evaluating the KT-equations ( ```kt_equations``` )

 [__NOTE__: Currently the only dispersion method is that in [[Dan14a]](https://arxiv.org/abs/1409.7708) but different methods such as a Pasquier-inversion method of [[Guo14]](https://arxiv.org/abs/1412.3970) can easily be implemented with a new namespace].

 #### Examples
The __kt_omega_unsubtracted.cpp__ script gives an example calculation of the _unsubtracted_ KT equations for the omega meson case for arbitrary iterations of the integral equation. Building is identical to above.

Because there is only one relevent isobar and no isospin dependence, we use the `isobar` class directly but for more complicated cases we will use a combination of `isobar` objects in the `kt_amplitude` class.

Subtracted dispersion relations will soon be fittable to data using the above `dalitz_fit` class.
