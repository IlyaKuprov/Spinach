# examples/dnp_sol/cavity_dnp/validate_cavity_dnp.m

Runs the native Spinach donor-cavity example with three, four, and five Fock states, confirming the declared Cn subsystem and full Hilbert dimensions 60/80/100. Checks cutoff convergence and density trace, Hermiticity, and positivity at absolute 1e-8 accuracy.

Independent native `evolution` controls compare a genuine half-step trajectory at every original sample time, zero coupling, lossless Rabi oscillations and excitation conservation, and exact damped Jaynes-Cummings survival and photon transients. The existing native generators are reused, rather than reconstructing a parallel dissipator.

Call `validation=validate_cavity_dnp()` from the example directory after configuring the Spinach path. The returned structure contains numerical errors and physical observables; failed tolerances raise assertions. A success marker is printed after the checks. No files are written, and the checks validate the stated single-donor idealisation rather than omitted ensemble physics.
