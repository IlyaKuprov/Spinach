# examples/spin_chemistry/singlet_yield_1.m

- Signature: `singlet_yield_1()`
- Source: [examples/spin_chemistry/singlet_yield_1.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/spin_chemistry/singlet_yield_1.m)
- Recombination callback: [experiments/spin_chem/rydmr_exp.m](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/spin_chem/rydmr_exp.m)

## Purpose

Computes a liquid-state magnetic-field effect on the singlet recombination yield of a radical pair containing two electrons and four protons. The example uses an exponential singlet-recombination model via `liquid(...,@rydmr_exp,...,'labframe')`; its callback source documents a singlet-singlet RYDMR model and constructs the initial two-electron singlet. The callback header points to [DOI 10.1080/00268979809483134](https://doi.org/10.1080/00268979809483134) as a model reference, not as evidence that this example's calculated values are experimental measurements.

## Spin system and sweep

- The isotope list is `{'E','E','1H','1H','1H','1H'}`; the electron Zeeman scalar entries are `2.002` for both electrons. The scalar-coupling matrix is passed as `mt2hz(matrix/2)`: the nonzero electron–proton entries are `0.195` for electron 1 with protons 3 and 4, `-1.3` for electron 2 with proton 5, and `0.2` for electron 2 with proton 6. The source gives no separate unit labels for those matrix entries before conversion.
- The primary magnet is set to `1` for normalisation. The callback documents field values in tesla and rate constants in Hz. The field grid is `1e-3*(0:0.01:5)` T (0–5 mT in 0.01 mT increments). Eight singlet-recombination rates are swept: `[0.176 0.880 1.76 3.52 8.8 17.6 35.2 52.8]*1e6` Hz.
- The callback evaluates the radical-pair singlet recombination model for electron indices `[1 2]`; the example requests the Zeeman operator and leaves ZTE off by default. Its basis is `sphten-liouv`, approximation `none`, zero-order projections `{0}`, with S2 symmetry on proton spins `[3 4]`.

## Observable and plot

The simulation returns `M`, plotted against the field grid with axes labelled magnetic field in tesla and singlet recombination yield. The source comment gives a calculation time of seconds. This is a computed field/rate sweep; the script does not report measured yield values.
