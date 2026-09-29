# examples/singlet_states/decoherence_bicyclopropylidene.m

[Source](https://github.com/IlyaKuprov/Spinach/blob/main/examples/singlet_states/decoherence_bicyclopropylidene.m)

This example constructs an eight-proton bicyclopropylidene spin system from vacuum-DFT coordinates, chemical shifts, couplings and chemical-shift-anisotropy tensors. The source describes a 65,536-dimensional Liouville space and includes every dipolar coupling and CSA tensor in the relaxation superoperator. Its `focus` is the spectrum of relaxation modes, rather than a simulated singlet-preparation or storage sequence.

The model selects 1H spins from `../standard_systems/bicyclopropylidene.log` (the conversion call also passes the numeric argument `31.8`, without a unit stated in the source), sets the field to 1.0 T, and uses Redfield relaxation, zero equilibrium, lab-frame retention and a 100 ps correlation time. The complete `sphten-liouv` basis is used without approximation; the relaxation-integration and zero tolerances are both 1e-5.

After constructing the relaxation superoperator, the function displays its 20 smallest-magnitude relaxation eigenvalues, labelled in Hz. It does not report an evaluated lifetime or a time-domain decay trace. The function contains no RF preparation pulse, gradient, storage interval or imaging reconstruction.
