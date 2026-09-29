# examples/spin_chemistry/singlet_yield_anisotropy_3.m

Source: [examples/spin_chemistry/singlet_yield_anisotropy_3.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/spin_chemistry/singlet_yield_anisotropy_3.m)

## Purpose

Compute the angular dependence of singlet yield for a model radical-pair reaction using Haberkorn recombination. The source comments estimate minutes on an NVIDIA Titan V and hours on a CPU; they are source estimates, not timings validated here. The line enabling GPU arithmetic is commented out.

## Spin system and interactions

The system is two electrons, two 14N nuclei, and three 1H nuclei. The field is set to 50e-6 T, labelled in the source as the Earth field. In the sphten-liouv basis, the source sets approximation IK-0 and inter_level=5. It assigns hyperfine tensors at (1,3), (2,4), (2,5), (2,6), and (2,7). The first two are passed to mt2hz as [-0.0989 0.0039 0; 0.0039 -0.0989 0; 0 0 1.7569] and [-0.0336 0.0924 -0.1354; 0.0924 0.3303 -0.5318; -0.1354 -0.5318 0.6680]. The remaining three use the same tensor, [-0.9920 -0.2091 -0.2003; -0.2091 -0.2631 0.2803; -0.2003 0.2803 -0.5398]. The source specifies scalar electron Zeeman values 2.0023 and 2.0025; nuclear values are zero. The matrices are source inputs to mt2hz; no additional unit label is supplied here.

## Recombination and angular yield

The reaction model is inter.chem.rp_theory='haberkorn' for electrons [1 2], with inter.chem.rp_rates=[1e6 1e6]. The powder calculation uses Lebedev grid leb_1ang_rank_63, tol=1e-2, verbose=0, and sum_up=0, and calls powder with @rydmr in the lab frame. The caller does not explicitly define the initial state; it delegates the singlet-yield computation to that routine.

The plotted output is cell2mat(yields) against grid.betas, labelled as the beta spherical angle in radians, with singlet yield on the vertical axis. The file identifies only a model reaction: it does not name the radicals or supply a measured yield.
