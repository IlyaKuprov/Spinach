# examples/esr_liq_pulsed/endor_phenyl.m

- Signature: `endor_phenyl()`

## Purpose

Simulate liquid-state Mims ENDOR of the phenyl radical. The isotropic g-factor (2.0024) and proton hyperfine couplings—ortho 17.4 G, meta 5.9 G, and para 1.9 G—are taken from Kasai, Hedaya, and Whipple (J. Am. Chem. Soc. 1969, 91, 4364). The two ortho protons and the two meta protons form equivalent pairs treated with S2 x S2 symmetry.

## Physical / mathematical content

- The spin system contains one electron and five ring protons: two ortho, two meta, and one para.
- The magnetic field is 0.33 T. Isotropic electron–proton hyperfine couplings are converted from 1.74, 0.59, and 0.19 mT to Hz using `mt2hz` and the phenyl-radical g-factor.

## Numerical / algorithmic content

- The basis uses `sphten-liouv` formalism with `approximation='none'`; `bas.sym_group` and `bas.sym_spins` specify the two equivalent proton pairs.
- Mims ENDOR is simulated with `liquid(spin_system,@endor_mims,parameters,'esr')`. Parameters specify zero offset, 512 points, a 300 MHz sweep, `tau=100 ns`, 4096-point zero filling, electron detection, and MHz axis units.
- The mean is subtracted from the simulated FID, which is then apodised with a Kaiser window of parameter 6. The code applies `fft`, centers the result with `fftshift`, and plots the spectrum magnitude against nuclear frequency in MHz.

## Implementation structure

- Define the field, isotopes, isotropic g-factor, and hyperfine couplings; configure the basis and symmetry; create the spin system; run the ENDOR simulation; process and plot the spectrum.
- The source comments give a calculation time of seconds.