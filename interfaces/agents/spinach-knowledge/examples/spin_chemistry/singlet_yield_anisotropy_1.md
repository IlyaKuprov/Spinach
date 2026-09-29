# examples/spin_chemistry/singlet_yield_anisotropy_1.m

Source: [examples/spin_chemistry/singlet_yield_anisotropy_1.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/spin_chemistry/singlet_yield_anisotropy_1.m)

## Purpose

Calculate the orientation-dependent singlet yield of a radical pair with the exponential recombination treatment selected by the source. The source comment estimates a run time of seconds; that estimate is not a measurement from this drafting pass.

## Spin system and interactions

The three-spin system is two electrons (E, E) and one proton (1H). The electron scalar Zeeman values are 2.0023 each and the proton value is 0. The only coupling tensor assigned is between spin 1 and spin 3; its diagonal input to gauss2mhz is [5e6 0 0; 0 4e6 0; 0 0 10e6]. The source does not state a unit for those input numbers, so they are reported as coded rather than relabelled. The basis is the full zeeman-hilb basis (bas.approximation='none').

## Powder calculation and observable

The caller sets sys.magnet=1 (commented as the unit-magnet field-sweep setting) and invokes powder with @rydmr_exp in the lab frame. Its sequence parameters are fields=50e-6, rates=2e6, electrons=[1 2], and Lebedev grid leb_2ang_rank_71; it requests zeeman_op, selects electron spins with spins={'E'}, and sets sum_up=0. The caller does not construct an initial density operator itself: that state and the detailed observable/kinetics handling are delegated to rydmr_exp.

The returned yield is converted from cells to an array and centred by subtracting its grid-weighted mean. The script maps that anisotropic yield over grid.betas and grid.gammas onto Cartesian coordinates, then colours a triangulated surface by the radius sqrt(x.^2+y.^2+z.^2). This is a visualisation of the calculated angular dependence, not a reported experimental yield.
