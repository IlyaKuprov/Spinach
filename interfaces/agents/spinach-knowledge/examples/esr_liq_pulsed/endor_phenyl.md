# examples/esr_liq_pulsed/endor_phenyl.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/esr_liq_pulsed/endor_phenyl.m`
- Signature: `endor_phenyl()`
- Total lines: 69

## Purpose

Mims ENDOR on a phenyl radical in liquid state. The g-factor and the isotropic proton hyperfine couplings are specified explicitly from the experimental data of Kasai, Hedaya, and Whipple (J. Am. Chem. Soc. 1969, 91, 4364): a(ortho)=17.4 G, a(meta)=5.9 G, a(para)=1.9 G, and an isotropic g-factor of 2.0024. The two ortho and the two meta protons are magnetically equivalent, so the full symmetry treatment uses an S2 x 

## Physical / mathematical content

- Liquid-state ESR examples. The dominant physics is electron Zeeman interaction, hyperfine coupling, relaxation broadening, and pulse-acquire or ENDOR-type detection in fast tumbling systems.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- Mims ENDOR on a phenyl radical in liquid state. The g-factor and the
- isotropic proton hyperfine couplings are specified explicitly from the
- experimental data of Kasai, Hedaya, and Whipple (J. Am. Chem. Soc.
- 1969, 91, 4364): a(ortho)=17.4 G, a(meta)=5.9 G, a(para)=1.9 G, and an
- isotropic g-factor of 2.0024. The two ortho and the two meta protons
- are magnetically equivalent, so the full symmetry treatment uses an
- S2 x S2 group direct product.
- Calculation time: seconds
- Magnet field
- Electron and the five ring protons (two ortho, two meta, one para)
- Phenyl radical isotropic g-factor
- Isotropic proton hyperfine couplings, converted from milliTesla

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `mt2hz()`, `create()`, `basis()`, `liquid()`, `apodisation()`, `fftshift()`, `kfigure()`, `plot_1d()`.
