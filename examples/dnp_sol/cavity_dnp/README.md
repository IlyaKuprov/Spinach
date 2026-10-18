# Minimal explicit-cavity DNP cycle

Time-domain illustration of one pumping cycle in Lorriaux et al., *Cavity-assisted dynamic nuclear polarisation*, arXiv:2610.19180v1. This is a single-donor model, not a reproduction or fit of the ensemble experiment.

After adding the Spinach root to the MATLAB path, run `fix_path('add')`, change into this directory, and run:

```matlab
r=cavity_dnp_minimal(3);
plot(r.time,r.populations);
xlabel('time / s'); ylabel('population');
legend('|6>','|7>','|15>');
validation=validate_cavity_dnp();
```

The example diagonalises the 20-level Bi:Si Hamiltonian with the paper's rounded parameters, retaining states |6>, |7>, and |15> in ascending-energy order. Seven spectator ground populations remain at 1/9 each; the propagated active density has trace 2/9 and the full density has trace one. All upper-manifold states start empty.

A selective microwave pi pulse swaps |6> and |15>, then a radiofrequency pi pulse swaps |6> and |7>. During the waiting stage the explicit truncated Weyl cavity couples |15,0> to |6,1>, and photon leakage returns population to |6,0>. The Liouville generator includes the photon-loss dissipator, not an eliminated Purcell-rate replacement. Coupling to the cavity is switched off during the ideal selective pulses; photon loss remains active.

The cavity width is kappa/(2*pi)=210 kHz. The one-donor coupling g/(2*pi)=100 Hz and pulse Rabi frequencies of 10 kHz are illustrative, not fitted experimental parameters. The waiting interval is eight nominal lifetimes using Gamma_P=4*g^2/kappa; this rate only sets the time window and the comparison curve. The near-twofold enhancement is the population ratio of |6> after one ideal cycle, not a collective enhancement or a twofold nuclear polarisation. The separate delta_iz output includes the physical nuclear expectation values of the three active eigenstates.

`validate_cavity_dnp` checks Fock cutoffs 2/3/4, trace, Hermiticity, positivity, the exact damped Jaynes-Cummings survival, half-step propagation, zero coupling, lossless Rabi oscillations, excitation conservation, and cavity-photon transients. It returns metrics without writing files. Both functions require Spinach and MATLAB; no parallel pool is required.

Excluded physics: collective emission, donor coupling/disorder distributions, experimental adiabatic pulse profiles, multiple cycles, spin diffusion, spectator dynamics, and non-cavity relaxation.
