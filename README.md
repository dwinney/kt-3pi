# kt_3pi
Codebase with tools to study J^PC -> 3 pi decays of mesons and rescattering effects with Khuri-Treiman formalism.

Requires ROOT (tested with version 6.17) with the ROOT::Minuit2 and MathMore libraries installed.
MathMore additionally requires a distribution of GSL (version >= 1.8) linked with ROOT.

## Dalitz Plot Parameter Extraction
The `dalitz_fit` object allows an arbitrary amplitude to be fit to the a polynomial expansion around the center of the Dalitz plot. __BW.cpp__ contains an example with a Breit-Wigner amplitude used to extract plot parameters for the omega meson decaying to three pions.   
To build in command-line use:
```
mkdir build && cd build
cmake ..
make BW && ./BW
```

Additionally you can make your own amplitude by inheriting and overriding virtual functions in the ```amplitude``` class and use the `dalitz` and `dalitz_fit` objects to extract parameters and make Dalitz plots. Changing the parameters stored in the ```decay_kinematics``` object allows you to consider arbitrary decay mass and quantum numbers.

## Khuri-Treiman Equations

The ```kt_amplitude``` class will allow to create dispersive isobar amplitudes within the KT formalism for the three pion final state.

 [__NOTE__: Currently implemented is the "standard" method of evaluating the dispersion relation up to infinity with a smoothly extrapolated phase and the conformal mapping scheme used in [[Dan14a]](https://arxiv.org/abs/1409.7708). This is changed with the `use_conformal` flag in `kt_options`, see examples below. Different methods such as a Pasquier-inversion method of [[Guo14]](https://arxiv.org/abs/1412.3970) can be later implemented].

 The `kt_options` class allows users to set and pass certain options to be used by the KT equation solver. These include:
 1. `int max_iter` sets the maximum number of times to iterate the KT solution.
 2. `int max_spin` sets the maximal spin to consider in the truncated isobar sum (currently numerical singularities appear in spin greater than 1, will be addressed in the future).
 3. `bool use_conformal` sets the method to evaluate dispersion relations. If `FALSE`, defaults to "standard" method of evaluating the dispersion relation up to infinity with a smoothly extrapolated phase and complex subtraction constants. If `TRUE` uses the conformal mapping scheme used in [[Dan14a]](https://arxiv.org/abs/1409.7708) with real subtraction constants.
 4. `double interp_cutoff` sets the maxminal invariant mass [GeV^2] which to save interpolations of the KT solution.


 #### Executables
To build any of the example executables simply use ``make`` in above.

__once_subtracted.cpp__ contains an example calculation of the KT equations for the omega meson case for arbitrary iterations of the integral equation to compare the effects of rescattering with only one subtraction constant which is fixed by fixing the total decay width to the experimental value.
