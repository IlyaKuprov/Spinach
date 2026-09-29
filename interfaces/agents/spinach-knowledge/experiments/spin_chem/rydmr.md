# experiments/spin_chem/rydmr.m

Source: [experiments/spin_chem/rydmr.m](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/spin_chem/rydmr.m)
Wiki: [rydmr.m](https://spindynamics.org/wiki/index.php?title=rydmr.m)

Signature: `A=rydmr(spin_system,parameters,H,R,K)`

This is a singlet-singlet radical-pair recombination calculation using the full chemical-kinetics superoperator. The inputs `H`, `R`, and `K` are same-sized numeric matrices: respectively the zero-external-field Hamiltonian commutation superoperator, relaxation superoperator, and chemical-kinetics superoperator. The source composes `L=H+1i*R+1i*K`.

The initial state is the two-electron singlet built from `spin_system.chem.rp_electrons(1:2)` and normalised to unit 2-norm. The solver computes the singlet projection with BICG and weights it by the first radical-pair recombination rate, using the source expression `A=rp_rates(1)*imag(S'*bicg(L,S,tol,numel(S)))`. The output `A` is a scalar fractional singlet yield.

The required `parameters` field is `tol`, a positive real scalar BICG tolerance. The source comment says `1e-2` is generally a good tolerance; this is guidance in the source, not a result measured for a particular system. The radical-pair electron indices and rate are taken from `spin_system.chem.rp_electrons` and `spin_system.chem.rp_rates`.
