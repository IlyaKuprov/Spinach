# experiments/traject.m

Source: [experiments/traject.m](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/traject.m)
Wiki: [traject.m](https://spindynamics.org/wiki/index.php?title=traject.m)

Signature: `traj=traject(spin_system,parameters,H,R,K)`

This function records forward time evolution of an initial state. It forms `L=H+1i*R+1i*K` from same-sized numeric matrix inputs `H`, `R`, and `K`. Its required parameter fields are `sweep`, `npoints`, `rho0`, and `decouple`. `sweep` is a positive real scalar sweep width in Hz; `npoints` is a positive integer; `rho0` is the initial state; and `decouple` is a cell array of isotope labels.

The time spacing is `1/sweep` seconds. The source calls `evolution` in `trajectory` mode for `npoints-1` propagation steps from `rho0` and describes the result as a system trajectory, a bookshelf stack of state vectors. Thus the requested trajectory has `npoints` states, beginning with the supplied initial state; it is not an acquired FID.

If `decouple` is non-empty, the requested isotope labels must occur in the spin system, and analytical decoupling is restricted to the `sphten-liouv` formalism. The routine applies the decoupling transformation to both `L` and `rho0` before propagation. The source example is `{'15N','13C'}`; an empty decoupling list skips that transformation.

This describes the source's propagation request and documented trajectory representation; no simulation was run to materialise a trajectory.
