# etc/gui/IO Formats/Gauss.cpp

[Implementation](https://github.com/IlyaKuprov/Spinach/blob/main/etc/gui/IO%20Formats/Gauss.cpp) · [Gaussian interface for MATLAB](https://github.com/IlyaKuprov/Spinach/blob/main/interfaces/gaussian/gparse.m)

The Windows reader builds one spin system rather than a calculation trajectory. It takes the last complete, finite input, standard or Z-matrix **orientation table** (Cartesian coordinates in ångströms). If no such table was printed it tries the Cartesian Gaussian archive, an explicitly labelled standalone thermochemistry Cartesian table, and the echoed atom list; internal-coordinate input is not a Cartesian geometry. Gaussian centre numbers are mapped explicitly: dummy atoms, ghost basis centres and translation vectors are not spins. Magnetic isotope labels take precedence over a geometry mass; if neither names an available isotope, the GUI uses the element's catalogue default. Every atomic number from 1 to 118 has a geometry entry; an unknown nuclear constant is not inferred from a site's coordinates. A new centre list in a linked calculation discards the old molecule's properties; a new spin multiplicity discards its previous EPR data.

| Printed quantity | GUI representation |
|---|---|
| Magnetic shielding | Nuclear tensor (all nine printed components are read; the GUI stores their symmetric part) in ppm, not a chemical shift. Use the component labels, not the order of the columns. For repeated tables use the last result, including final composite ONIOM shielding. |
| Total nuclear spin–spin J | One isotropic spin pair in Hz per off-diagonal entry; rows span successive column panels on larger systems. |
| Electronic g | Dimensionless electron tensor. For a hyperfine-only open-shell system the GUI supplies the free-electron g value. |
| Hyperfine | Electron–nucleus tensor in gauss. Add the final isotropic contact to the lab-frame dipolar tensor; reconstruct that tensor from Baa/Bbb/Bcc, their printed direction cosines and `multiplicity − 1`. The centre and isotope are on the Bbb row. In ONIOM jobs use the integrated real-system contact and anisotropic spin densities, scaled by the nuclear magnetic moment and rotated into the stored coordinate frame. Link atoms from model-layer calculations do not describe the corresponding real nuclei. |
| Nuclear quadrupole coupling [Chi] | Traceless nuclear tensor in MHz. Its components are in the molecular principal-axis frame; align the corresponding coordinates before storing it in the lab frame. For a linear molecule the printed principal basis distinguishes the two transverse axes. |
| Molecular susceptibility | Separate molecular tensor in the printed cgs-ppm convention. It is not an atom, electron or interaction-grid row; inspect the susceptibility display or Spinach export. |

A missing property section does not invalidate a geometry. Gaussian overflow markers are not physical zeros. Orbital energies, vibrational spectra, density matrices and optimisation histories are not imported.
