# Native Spinach cavity-assisted DNP cycle

Time-domain illustration of one ideal pumping cycle in Lorriaux et al., *Cavity-assisted dynamic nuclear polarisation*, arXiv:2610.19180v1. This is a single-donor model, not an ensemble/collective calculation or experimental fit.

After adding the Spinach root to the MATLAB path, call `fix_path('add')`, change into this directory, and run:

```matlab
r=cavity_dnp_minimal(3);
plot(r.time,r.populations);
xlabel('time / s'); ylabel('population');
legend('|6>','|7>','|15>');
validation=validate_cavity_dnp();
```

## Native infrastructure

`sys.isotopes={'E','209Bi',['C' num2str(nlevels)]}` declares the electron, bismuth nucleus, and explicit cavity. `create`, `basis`, and `hamiltonian` construct the full donor-cavity system. All 20 donor levels and every retained Fock level remain in the density matrix: its Hilbert dimension is `20*nlevels`, not three donor states with an externally appended cavity. Spinach requires at least three levels for a native cavity.

`operator` supplies donor observables, cavity annihilation/number operators, and Fock projectors; `state` supplies the cavity vacuum. The donor Hamiltonian vacuum block is diagonalised to identify the paper's ascending-energy state labels. Spectral projectors select the two ideal transitions and are extended to occupied cavity levels using the native creation operator. This spectral selection specifies the experiment; it is not a separate state-space or propagator implementation.

`hilb2liouv` constructs the selected-transition commutators. `inter.modes.linewidths` and `relaxation` construct the native cavity GKSL loss. `evolution` propagates both pulses and the emission trajectory. No hand-written Kronecker commutator, dissipator, or time-stepping loop is used.

## Physics and limits

The field is 65.37 mT and the isotropic hyperfine coupling is 1.475 GHz. Spinach supplies the electron and 209Bi magnetogyric ratios. Nine lower-manifold levels initially have population 1/9; all upper levels and non-vacuum Fock states start empty. The full density has trace one, including the seven spectator ground levels.

A selective microwave pi pulse swaps |6> and |15>, then a radiofrequency pi pulse swaps |6> and |7>. During the waiting stage the resonant selected-transition Jaynes-Cummings Hamiltonian couples |15,0> to |6,1>, and native photon loss transfers population to |6,0>. The donor-cavity coupling is switched off during the ideal pulses, while photon loss remains active. The selected resonant Hamiltonian is expressed in its interaction frame; the GHz donor Hamiltonian is used for physical state identification, not propagated in that frame.

The cavity width is kappa/(2*pi)=210 kHz. The one-donor coupling g/(2*pi)=100 Hz and pulse Rabi frequencies of 10 kHz are illustrative. The waiting interval is eight nominal lifetimes, Gamma_P=4*g^2/kappa; this rate sets the time window only, and is not substituted for the explicit cavity dynamics. Near-twofold enhancement means the final |6> population divided by 1/9, not a collective enhancement or doubled nuclear polarisation. `delta_iz` reports the physical nuclear expectation change separately.

Excluded physics: collective emission, donor coupling/disorder distributions, experimental adiabatic pulse profiles, repeated pumping, spin diffusion, and non-cavity relaxation.

## Validation

`validate_cavity_dnp` runs native Fock dimensions 3/4/5 (full Hilbert dimensions 60/80/100). It checks trace, Hermiticity, positivity to absolute 1e-8 numerical accuracy, cutoff convergence, exact damped Jaynes-Cummings survival, a genuine half-step trajectory at all original sample times, zero-coupling behaviour, lossless Rabi oscillations and excitation conservation, and resolved cavity-photon transients. The retained full-system generators allow these independent controls without rebuilding any dissipator. No files or figures are written by the two functions.
