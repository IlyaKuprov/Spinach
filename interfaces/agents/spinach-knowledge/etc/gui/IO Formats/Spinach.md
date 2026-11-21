# etc/gui/IO Formats/Spinach.cpp

[Implementation](https://github.com/IlyaKuprov/Spinach/blob/main/etc/gui/IO%20Formats/Spinach.cpp)

The exporter writes a MATLAB function named after the chosen file, not a MAT-file or a spectrum. It constructs sys and inter, then emits create(sys,inter) and basis(spin_system,bas). Nuclear coordinates are rows in inter.coordinates; electronic entries are empty. Chemical shifts and g tensors go to inter.zeeman.matrix. Scalar J goes to inter.coupling.scalar; hyperfine, exchange, ZFS and quadrupolar matrices go to inter.coupling.matrix. Display bonds are not pair Hamiltonian terms.

The writer rotates interaction matrices to the lab frame. It removes one third of the trace from each quadrupolar diagonal before formatting, converts frequency-valued couplings to Hz, and emits a gauss2mhz conversion for Gauss hyperfine data. Molecular susceptibility, when present, is written separately as inter.suscept.chi in cubic angstroms and inter.suscept.xyz at [0 0 0]. Its scale factor is 4π·10¹⁸/N_A for input in cgs-ppm. Check a nonzero off-diagonal and the centre as well as the units.
