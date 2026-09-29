# kernel/assume.m

Source: [kernel/assume.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/assume.m)

Signature: `spin_system = assume(spin_system,assumptions,retention)`

Selects interaction-strength rules for a simulation context. Call it before requesting the Hamiltonian. It writes the chosen assumptions string and strength selectors into spin_system.inter; it does not construct a Hamiltonian, a state vector, or Hilbert-space operator matrices.

## Selector layout

For N = spin_system.comp.nspins, the function replaces spin_system.inter.zeeman.strength and spin_system.inter.giant.strength with N×1 cell arrays, and spin_system.inter.coupling.strength with an N×N cell array. The cells contain rule strings such as secular, zz, full, strong, z*, *z, or ignore. These are selectors for later Hamiltonian generation, not operators. There are no state/operator dimensions or physical units defined by this function. It does not read or write a cache.

## Assumption sets

- nmr: high-field NMR rotating-frame rules. Zeeman and giant-spin selectors are secular; same-isotope couplings are secular and different-isotope couplings use zz. This set rejects component types C, V, or T.
- esr and deer: the same branch, with electrons in the rotating frame and nuclei in the laboratory frame. Electron-nucleus couplings use the directional z* or *z selector, electron-electron couplings are secular, and nuclear-nuclear couplings are strong.
- deer-zz: DEER selectors with electron-electron flip-flop terms removed using the zz coupling selector.
- labframe: full laboratory-frame spin interactions; bosonic modes are allowed and their mode terms remain in the laboratory frame.
- qnmr: spin-1/2 nuclei use rotating-frame rules and higher-spin nuclei use laboratory-frame rules. It rejects component types C, V, or T.
- se_dnp_h+, se_dnp_h-, and se_dnp_h0: select the corresponding positive-, negative-, or zero-frequency solid-effect DNP component. The source retains the matching EzNp, EzNm, or EzNz electron-nuclear terms and the specified T(L,+1), T(L,-1), or secular inter-nuclear terms, while ignoring inter-electron, giant-spin, quadratic, and Zeeman terms. These sets reject component types C, V, or T.
- cavity: a common rotating frame with the RWA; spin Zeeman, giant-spin, and spin-spin selectors are secular. When modes are present, the source sets mode frequencies to offsets, anharmonicity and Kerr terms to full, exchange to rwa, dispersive terms to full, and longitudinal/modulation terms to ignore. It checks that the longitudinal and modulation input terms being dropped are absent.
- spin-phonon: electrons use a rotating frame, while nuclei and bosonic modes remain in the laboratory frame. Electron Zeeman selectors are secular and nuclear Zeeman selectors full; electron-nucleus couplings use z* or *z. When modes are present, the source sets exchange to nonelec and Kerr, longitudinal, dispersive, coupling-modulation, and Zeeman-modulation selectors to full.

The assumptions argument must be a character string, and the source accepts two or three inputs; an unrecognised assumption string errors; the optional retention string can be couplings (sets all Zeeman selectors to ignore) or zeeman (sets all giant-spin and pair-coupling selectors to ignore). The retention pass changes spin selectors only; the source documents these filters as undefined for systems containing bosonic modes.

## Reference

[Spin Dynamics Wiki: assume.m](https://spindynamics.org/wiki/index.php?title=assume.m)
