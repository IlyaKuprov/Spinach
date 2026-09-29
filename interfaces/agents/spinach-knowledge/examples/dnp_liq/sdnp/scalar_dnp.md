# examples/dnp_liq/sdnp/scalar_dnp.m

- MATLAB implementation: [examples/dnp_liq/sdnp/scalar_dnp.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/dnp_liq/sdnp/scalar_dnp.m)

- Signature: `scalar_dnp()`
- Calculation time: seconds
- [MATLAB source](../../../../../../examples/dnp_liq/sdnp/scalar_dnp.m)
- Simulation reference: [Journal of Magnetic Resonance Open (2022)](https://doi.org/10.1016/j.jmro.2022.100040)
- Experimental data: [Angewandte Chemie International Edition](https://doi.org/10.1002/anie.201811892)

## Purpose

Calculates the field dependence of the Overhauser coupling factor for the `13C` nucleus of CHCl3 coupled to a nitroxide electron, and overlays four experimental values. The script obtains nuclear longitudinal relaxation and electron-to-nuclear cross-relaxation projections from the relaxation superoperator, then plots their ratio `Rx/R1n`.

## Spin system and relaxation settings

The spins are `E` and `13C`, with dipolar coordinates (0, 0, 0) and (0, 0, 3.1) Angstrom. Their Zeeman tensors are electron [2.0029, 2.0065, 2.0098] as a dimensionless g tensor (the Bohr magneton enters separately in the spin Hamiltonian) and nuclear [100, 120, 120] ppm; both Euler-angle triples are zero radians. The static isotropic hyperfine coupling is 2e6 Hz. The basis is `sphten-liouv` with no approximation.

Three relaxation contributions are enabled: `SRFK`, `redfield`, and `t1_t2`. The empirical electron rates are R1=2e6 Hz (T1=500 ns) and R2=5e6 Hz (T2=200 ns); the nuclear rates are R1=0.25 Hz (T1=4 s) and R2=3.00 Hz (T2=1/3 s). Rotational Redfield relaxation uses a 30 ps correlation time. Scalar collisional Redfield uses weighted correlation-time components (0.62, 30 ps) and (0.38, 0.80 ps), with scalar modulation depth 3.6e6 Hz for the pair. Temperature is 298 and equilibrium is Di Bari; relaxation retention is secular.

## Field sweep and rate ratio

The field values labelled in Tesla are `[logspace(-2,1,30),15,20,25,30]`: 30 logarithmic points from 0.01 to 10, followed by 15, 20, 25, and 30. A `parfor` loop creates the spin system at each field and obtains `R=relaxation(spin_system)`. The script normalises the `13C` and electron `Lz` states, then calculates `R1n=real(Nz'*R*Nz)` and `Rx=real(Nz'*R*Ez)`; the plotted coupling factor is `Rx./R1n`.

The four tabulated experimental pairs (field in Tesla, coupling factor) are (0.34, −0.84), (1.20, −0.48), (3.40, −0.37), and (9.40, −0.16). The plot uses a logarithmic field axis, shows the calculated curve and data markers together, and sets the displayed axes to 0.01–30 Tesla and −1 to 0. This is a relaxation-rate projection for the specified two-spin model, not a propagated microwave-driven trajectory.
