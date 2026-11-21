# etc/gui/IO Formats/Castep.cpp

[Implementation](https://github.com/IlyaKuprov/Spinach/blob/main/etc/gui/IO%20Formats/Castep.cpp)

The CASTEP reader accepts modern MAGRES v1 atom and magnetic-property blocks (square or angle brackets), legacy MAGRES text tables and printed CASTEP geometries. A site is identified by its label and numeric index, not its position in a property table. Joined identifiers such as O100 are matched against the atom records. Delimited identifiers remain distinct even when their joined spellings collide. The order of geometry and magnetic-property blocks does not determine which atom receives a tensor.

Cartesian MAGRES coordinates are in angstroms. In a .castep report the reader takes the last complete fractional-coordinate table and converts it using that frame's real lattice vectors. A report that only gives shielding invariants, without the full oriented tensor, is not equivalent to a .magres file.

Total ms is chemical shielding in ppm; the GUI tensor model stores the symmetric part of the printed matrix. Total efg is in atomic units, not MHz. The GUI converts it to a quadrupolar matrix in MHz using the assigned isotope's spin and quadrupole moment. If necessary it selects a quadrupolar isotope of the same element from its isotope table. There is no quadrupole interaction for spin below one. Legacy tables supply their own quadrupole moments.

Total isc is a reduced coupling, not J in Hz. The reader averages reciprocal records' isotropic parts and converts to Hz with the signed gyromagnetic ratios of the two assigned isotopes. Self records and component decompositions are not additional couplings. Energies, bands, forces and the periodic sample's macroscopic susceptibility are not molecular spin-system properties.
