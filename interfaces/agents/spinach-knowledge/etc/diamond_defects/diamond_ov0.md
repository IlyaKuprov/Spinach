# etc/diamond_defects/diamond_ov0.m

- MATLAB implementation: [etc/diamond_defects/diamond_ov0.m](https://github.com/IlyaKuprov/Spinach/blob/main/etc/diamond_defects/diamond_ov0.m)

## Purpose and call

`[sys,inter]=diamond_ov0(parameters)` constructs the ground-state spin system for the neutral oxygen-vacancy centre OV0 (WAR5) in diamond. The input is a structure with a character-string `orientation` of `'111'`, `'110'`, or `'100'`; the selected crystal-plane normal is aligned with the magnetic field. There is no default orientation, and this function does not take a field magnitude.

## Physical model and parameters

The defect has electron spin S=1 and C3v symmetry. Its electron g and zero-field-splitting tensors are axial about the trigonal axis; rhombicity is zero within experimental error. The g principal values are 2.0025, 2.0025, and 2.0029, and the zero-field splitting is 2888 MHz. Oxygen is omitted because predominant natural-isotope 16O has no nuclear spin.

The model covers four resolved 13C hyperfine shells, labelled a, g, l, and c, with 3, 6, 3, and 3 symmetry-equivalent atoms, respectively. It returns one representative 13C per shell (four carbons total), not all 15 sites. The other site tensors can be obtained by C3v operations about the trigonal axis. The a and g parameters use the symmetry-relaxed fit, reported in the thesis as the better fit.

| Shell | Hyperfine principal values (MHz) | Principal-axis polar angles from [001] (degrees) | Azimuths from [100] (degrees) |
|---|---|---|---|
| a | 197.4, 117.3, 118.2 | 53.86, 143.86, 90.00 | 225.26, 225.00, 135.00 |
| g | 17.5, 11.7, 13.0 | 60.40, 150.40, 90.00 | 225.26, 225.00, 135.00 |
| l | 12.6, 8.5, 8.5 | 54.70, 144.74, 90.00 | 225.26, 225.00, 135.00 |
| c | 7.4, 4.3, 4.3 | 54.70, 144.74, 90.00 | 225.26, 225.00, 135.00 |

The code forms and symmetrises hyperfine tensors from the principal values and crystal-frame axes, converts the MHz values to Hz for the interaction matrices, then rotates the tensors for the requested orientation.

## Outputs and limitations

`sys` contains isotope entries `E3` and four `13C` spins, labelled `OV0`, `C_a`, `C_g`, `C_l`, and `C_c`. `inter` contains the electron Zeeman tensor, electron zero-field-splitting tensor, and four electron–13C hyperfine tensors. The function returns specifications rather than a simulated spectrum; it does not explicitly include 16O or every symmetry-equivalent carbon site.

## Sources

- B. L. Cann, *Magnetic Resonance Studies of Point Defects in Diamond*, PhD thesis, University of Warwick (2009), Tables 9-2 and 9-3: <https://wrap.warwick.ac.uk/id/eprint/3125>
- S. Mukherjee et al., *Phys. Rev. B* **114**, 074105 (2026): OV0 assignment and zero-field-splitting confirmation at 4 K, <https://doi.org/10.1103/3dcd-mkcq>
- Source documentation: <https://spindynamics.org/wiki/index.php?title=diamond_ov0.m>
