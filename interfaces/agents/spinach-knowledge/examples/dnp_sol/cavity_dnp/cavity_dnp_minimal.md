# cavity_dnp_minimal

`result=cavity_dnp_minimal(nlevels)` models an ideal shelving, conditional nuclear swap, and cavity-mediated reset cycle using an electron, a spin-1/2 nucleus, and a Cn bosonic mode. It is a four-level high-field analogue of Fig. 1 in Lorriaux et al., arXiv:2610.19180v1, not a low-field Bi:Si donor simulation or an ensemble fit.

The formalism is exclusively `sphten-liouv`. Native interaction specifications generate the Hamiltonian and cavity relaxation; native state/operator syntax generates vacuum, spin populations, conditional pulse generators, and observables. `step` applies ideal pulses; `evolution` propagates explicit spin-cavity exchange and photon loss. Physical expectation values use native left-product operators contracted with the identity trace functional. Returned fields include populations, Iz, photons, trace, Fock-edge occupation, enhancement, and native generators for controls.

The initial nuclear mixture is unpolarized, the cavity bath is at zero temperature, and ideal pulses omit drift and decay. Illustrative rates are not the paper's experimental parameters. See the adjacent README for approximations and `validate_cavity_dnp` for numerical controls.
