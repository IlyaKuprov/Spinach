# examples/quantum_tech/spin_phonon_dephasing.m

- Signature: `spin_phonon_dephasing()`
- Source: [MATLAB example](https://github.com/IlyaKuprov/Spinach/blob/main/examples/quantum_tech/spin_phonon_dephasing.m)

## Model

This source presents longitudinal spin–phonon coupling as a minimal Weyl-algebra model for strain-modulated spin dynamics discussed in the context of NV centres and molecular spins. It does not instantiate a specific NV defect or molecular electronic structure. The isotope list `{'E','V7'}` means a multiplicity-2 electron (spin one-half) and a seven-level phonon mode; `V7` is the phonon truncation, not a defect or vanadium isotope. The source sets `sys.magnet=0`, phonon frequency `20e6` (20 MHz), and longitudinal coupling `4e6*sqrt(2)`. It specifies no EPR field/frequency selection, drive, relaxation, or dissipation.

## State, propagation, and observables

The calculation uses the `zeeman-hilb` formalism without a basis approximation. Its initial state is a sum of a transverse-spin `Lx` component and a `ZL2` population component, each paired with the phonon vacuum `BL1`. A single `spin-phonon` device trajectory is used. The source notes that the drift commutes with spin projection, so these coherence and population sectors evolve without mixing in this model.

At each trajectory point it projects the spin coherence `Lx` and the phonon displacement quadrature `C + A` (annihilation plus creation operators). The two real expectation traces are plotted over 0–1 μs with 501 points. The source includes checks that coherence modulation and oscillator displacement are visible; those checks do not establish a measured dephasing rate or an experimental defect spectrum.

## Interpretation

The plots illustrate modelled coherence modulation and spin-conditioned oscillator displacement under longitudinal coupling. The example supplies no measured parameters, stochastic bath, or defect-specific isotope selection, so the curves should not be read as measured NV-centre data or a prediction of device fidelity. Its “Calculation time: seconds” comment is not a runtime validation or convergence claim.
