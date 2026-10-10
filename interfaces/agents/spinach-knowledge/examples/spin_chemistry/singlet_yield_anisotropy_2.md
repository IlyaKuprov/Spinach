# examples/spin_chemistry/singlet_yield_anisotropy_2.m

Source: [examples/spin_chemistry/singlet_yield_anisotropy_2.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/spin_chemistry/singlet_yield_anisotropy_2.m)

## Purpose

Calculate a radical-pair singlet-yield anisotropy with exponential recombination kinetics. The source comment estimates a run time of seconds; it is not a measured timing from this drafting pass.

## Spin system and anisotropic couplings

The model has two electrons, two 14N nuclei, and one 1H nucleus, in that spin order. The scalar Zeeman values are 2.0023 for each electron and 0 for the three nuclei. Three anisotropic tensors connect spin pairs (1,3), (2,4), and (1,5). Their diagonal eigenvalue matrices are A1=diag(-1.049,-0.996,13.826), A2=diag(-0.305,-0.222,6.872), and A3=diag(-13.850,-9.372,0.143). The rotation matrices are R1=[0.4380 0.8655 -0.2432; 0.8981 -0.4097 0.1595; -0.0384 0.2883 0.9568], R2=[0.9703 -0.2207 0.0992; 0.2383 0.9426 -0.2340; -0.0419 0.2506 0.9672], and R3=[0.9819 0.1883 -0.0203; -0.0348 0.2850 0.9579; -0.1861 0.9398 -0.2864]. The source passes each rotated tensor through 1e6*gauss2mhz(R*A*R') into the coupling entries. These are the coded inputs and conversion; the source does not provide a unit label for the eigenvalue numbers.

The model uses the full zeeman-hilb basis (bas.approximation='none'). No molecular identity or explicit initial density operator is specified in this caller.

## Powder calculation and visualisation

The field-sweep setup uses sys.magnet=1, with sequence values fields=50e-6 and rates=50e6, electron indices [1 2], and Lebedev grid leb_2ang_rank_35. The powder calculation calls @rydmr_exp in the lab frame with spins={'E'}, needs={'zeeman_op'}, and sum_up=0; exponential kinetics and the yield calculation are handled through that routine rather than a separately constructed state in this script.

The script converts the cell output to an array and subtracts its grid-weighted mean. It places the centred yield on the spherical directions grid.betas and grid.gammas, then renders a triangulated, colour-mapped surface with trisurf and edge transparency 0.25. The surface encodes the computed anisotropy; no numerical or experimental yield is asserted here.
