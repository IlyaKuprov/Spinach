# examples/nmr_solids/fitting/bromide_csa_nqi/kbr_mas_fitting.m

- Signature: `kbr_mas_fitting()`

## Purpose

Fitting of a 79Br MAS NMR spectrum of potassium bromide with respect to the quadrupole coupling constant. The spectrum cannot be fitted with a single quadrupolar tensor; at least 3 are necessary, likely due to a dist- ribution of electrostatic environments in the powder. Calculation time: hours.

## Physical / mathematical content

- Solid-state NMR examples. The key physics is anisotropic spin interactions under static or magic-angle-spinning conditions: chemical-shift anisotropy, dipolar coupling, quadrupolar coupling, cross-polarisation, and orientation averaging using Floquet, Fokker-Planck, or direct powder quadrature formalisms.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.
- Chemical-shift anisotropy is present: shielding is treated as a second-rank tensor whose orientation relative to the field or rotor axis modulates line shapes and transfer dynamics.
- Quadrupolar physics is relevant: nuclei with spin > 1/2 interact with the electric field gradient tensor, introducing second-rank anisotropy, asymmetry, and overtone or MQ phenomena.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- Fitting of a 79Br MAS NMR spectrum of potassium bromide
- with respect to the quadrupole coupling constant.
- The spectrum cannot be fitted with a single quadrupolar
- tensor; at least 3 are necessary, likely due to a dist-
- ribution of electrostatic environments in the powder.
- Calculation time: hours.
- Load and normalise the data
- Set instrumental variables
- Set optimizer options
- Get a figure going
- Run the optimisation
- Plot and print the fitted parameters
- Least squares error function
