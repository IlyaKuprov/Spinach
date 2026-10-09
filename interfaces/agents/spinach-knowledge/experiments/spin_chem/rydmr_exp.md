# experiments/spin_chem/rydmr_exp.m

Source: [experiments/spin_chem/rydmr_exp.m](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/spin_chem/rydmr_exp.m)
Wiki: [rydmr_exp.m](https://spindynamics.org/wiki/index.php?title=rydmr_exp.m)
Reference: [DOI 10.1080/00268979809483134](https://doi.org/10.1080/00268979809483134).

Signature: `answer=rydmr_exp(spin_system,parameters,H,R,K)`

This routine computes singlet-singlet radical-pair yields over field and singlet-recombination-rate inputs with the exponential-recombination treatment described by the cited paper. It expects a unit primary magnet specification (`sys.magnet=1`), and the supplied system magnet must equal one. Its source note says this function supplies exponential recombination itself, so do not combine it with additional recombination loss records in `inter.chem.reactions`.

The parameter row vector `fields` contains magnetic fields in tesla; `rates` contains singlet recombination rates in Hz. `electrons` gives the two electron indices in the isotope list, with `[1 2]` as the source example. The function checks that they are two positive integer indices within the system and that both selected spins are electrons. The caller should request `'zeeman_op'` in `parameters.needs` so the context supplies the field-sweep operator `parameters.hzeeman`. The routine also takes same-sized numeric matrix inputs `H`, `R`, and `K`.

The function sets `Z=parameters.hzeeman` and removes that term from the supplied Hamiltonian before applying each field. For `sphten-liouv` and `zeeman-liouv` it forms `L=H+1i*R+1i*K`, then evaluates each rate-field pair with the field term and the rate-dependent exponential loss included in that pair's Liouvillian. For `zeeman-hilb` it requires `R` to be a scalar multiple of the identity; it separates the singlet's unit-state component from its traceless part, propagates the latter to the source-defined endpoint `t_end=10/rate`, and adds the unit-state contribution analytically. These are distinct source branches; unsupported formalism values produce an error.

The returned array is `reshape(real(answer),N,M)` with `N=numel(rates)` and `M=numel(fields)`: rows correspond to rate values and columns to field values. The source documents `fields` and `rates` as row vectors. The DOI, Wiki page, normalisation convention, electron-index example, and warning about separate chemical rates are retained from the source documentation.
