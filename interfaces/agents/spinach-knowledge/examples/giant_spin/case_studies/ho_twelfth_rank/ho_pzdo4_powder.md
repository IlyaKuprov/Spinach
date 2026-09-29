# examples/giant_spin/case_studies/ho_twelfth_rank/ho_pzdo4_powder.m

- MATLAB implementation: [examples/giant_spin/case_studies/ho_twelfth_rank/ho_pzdo4_powder.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/giant_spin/case_studies/ho_twelfth_rank/ho_pzdo4_powder.m)

- Signature: `ho_pzdo4_powder()`

## Purpose

Calculates powder-averaged pulsed-field magnetisation for the Ho(pzdo)4 metal-organic framework, modelled as a J=8 giant spin with 90 crystal-field coefficient slots through rank 12 (from the even-rank 2–12 parameter set). The setup follows the CASSCF parameter set and Fig. 3 case study in [arXiv:2609.16352](https://arxiv.org/abs/2609.16352).

## Physical model and setup

The Spinach system has one `E17` effective spin, g = 1.24, and crystal-field coefficients from `ho_pzdo4_params.m`, converted rank by rank from cm^-1 to Hz and then to spherical form; the crystal-field Euler angles are zero. The Zeeman-Hilbert-space formalism is used without approximation, at 2 K. The field sweep is linear from 0 to 10 T at 10 T/ms, represented as B(t) = 1e4 t with t in seconds. The simulation uses 100,000 steps of 10 ns, spanning 1 ms, and records 1,000 output intervals.

Spin-phonon relaxation uses a generalised Lindblad dissipator in the form of Saito and Miyashita. The coupling matrix has unit elements between adjacent mJ states, with projections rounded to half-integers. The bath is super-Ohmic (alpha = 2), with lambda = 10 cm^-1 and I0 = 1e-14 ps/rad; the source converts the bath parameter to rad/s units for Spinach. The magnetisation observable is the moment -1.24 Jz, in Bohr magnetons.

## Powder calculation and output

The source ramps the field with `field_prof(t)=1e4*t` to 10 T over 1 ms; it uses 10 ns steps (`timestep=1e-8`), `nsteps=1e5`, and `nout=1000`, recording every 1000 steps.

The `powder` calculation uses the two-angle Lebedev grid `leb_2ang_rank_23` (194 orientations, the nearest grid available here to the paper's order-21, 170-point grid). Each orientation is retained separately (`sum_up=false`): the plot shows its curve in grey, their weighted powder average in red, and the thermal-equilibrium powder average in blue. The equilibrium values are calculated at the same field points using the full and Zeeman Hamiltonians for each grid orientation.

The plotted axes span 0-10 T and 0-8 Bohr magnetons. The source configures 16 processes for the powder calculation. The plotting configuration follows the paper's Fig. 3 case study.
