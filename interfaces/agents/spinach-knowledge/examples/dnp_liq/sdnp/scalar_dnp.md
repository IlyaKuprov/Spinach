# examples/dnp_liq/sdnp/scalar_dnp.m

- Signature: `scalar_dnp()`
- Simulation reference: [*Journal of Magnetic Resonance Open* (2022)](https://doi.org/10.1016/j.jmro.2022.100040)
- Experimental data: [*Angewandte Chemie International Edition*](https://doi.org/10.1002/anie.201811892)
- Calculation time: seconds

## Purpose

Calculates the field dependence of the Overhauser coupling factor for a 13C nucleus and nitroxide electron, and compares the simulated curve with four experimental values. The observable is the cross-relaxation rate divided by the 13C longitudinal relaxation rate.

## Spin system and relaxation model

The two spins are an electron and 13C separated by 3.1 Å, with a 2 MHz isotropic hyperfine coupling. The source specifies anisotropic electron g and 13C chemical-shift tensors, Redfield and scalar-collision (SRFK) relaxation, and empirical T1/T2 rates. It uses a 30 ps rotational correlation time, a two-component SRFK correlation model with weights 0.62 and 0.38, and a 3.6 MHz scalar modulation depth. The temperature is 298 K, the equilibrium model is Di Bari, and the relaxation retention is secular. The basis is sphten-liouv without an approximation.

## Field scan and output

At each of 34 magnetic fields (0.01–10 T logarithmically, followed by 15, 20, 25, and 30 T), the script builds the system and relaxation superoperator. It normalises the 13C and electron longitudinal states, evaluates the 13C relaxation rate and electron-to-13C cross-relaxation rate, and forms their ratio. A parallel loop computes the field points. The plot uses a logarithmic field axis and overlays the four tabulated experimental points.
