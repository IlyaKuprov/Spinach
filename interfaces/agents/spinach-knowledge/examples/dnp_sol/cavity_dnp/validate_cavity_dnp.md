# examples/dnp_sol/cavity_dnp/validate_cavity_dnp.m

Runs the explicit-cavity DNP example with two, three, and four Fock states and checks density-matrix invariants and cutoff convergence. Independent controls test zero coupling, lossless Rabi oscillations and excitation conservation, exact damped Jaynes-Cummings survival and photon transients, and a halved propagation step.

Call `validation=validate_cavity_dnp()` from the example directory after configuring the Spinach path. The returned structure contains numerical errors and physical observables; failed tolerances raise an assertion. The function prints a success marker after all checks and writes no files. These checks validate the stated single-donor idealisation, not the omitted ensemble physics.
