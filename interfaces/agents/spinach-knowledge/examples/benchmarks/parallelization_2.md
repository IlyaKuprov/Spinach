# examples/benchmarks/parallelization_2.m

- MATLAB implementation: [examples/benchmarks/parallelization_2.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/benchmarks/parallelization_2.m)

`parallelization_2()` has no arguments. It constructs an eleven-spin pyrene-cation model from the vacuum-DFT log at the relative path ../standard_systems/pyrene_cation.log; gparse and g2spinach must be available, and that file must resolve from the run directory. The parser maps all ten hydrogens to `1H` and appends one `E` electron spin in EPR mode, with `options.no_xyz=1`.

The benchmark sets the field to 50 µT, uses the Zeeman Hilbert-space formalism with no basis approximation, assumes the lab frame, and adds the Hamiltonian contribution returned by orientation(Q,[pi/3,pi/4,pi/5]). Its initial operator is Lz on the electron spins. It calls evolution(...,5e-9,200,'observable') to time a 200-step observable propagation. The source does not label the 5e-9 argument's unit.

MATLAB's Parallel Computing Toolbox is required: candidate pool sizes are [2 4 8 16 32 64 128 256 512 1024], filtered to feature('numcores'). For each size the script deletes any current pool, starts a new parpool, waits 10 seconds, then times propagation and prints the elapsed seconds. Thus pool creation and the settling pause are outside the timed interval; the final pool remains active when the function exits.

The only reported output is one timing line per usable pool size; the script saves no data and makes no figure. Its DOI is [10.1063/1.3679656](https://doi.org/10.1063/1.3679656). Timings are machine- and MATLAB-dependent; the example itself does not establish a particular speedup or validate the propagated observable against a reference.
