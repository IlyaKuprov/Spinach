# etc/gui/IO Formats/Orca.cpp

[Implementation](https://github.com/IlyaKuprov/Spinach/blob/main/etc/gui/IO%20Formats/Orca.cpp)

ORCA Cartesian coordinates marked A.U. are scaled by 0.529177249 to angstroms. The g-matrix, zero-field splitting, chemical shifts and electric/magnetic hyperfine structure use separate sections. The importer represents the electronic centre explicitly; a property table need not have one row per visible nucleus.

Check the interaction's stored unit rather than treating every printed matrix as ppm. The shift and hyperfine paths use MHz; ZFS uses inverse centimetres. The quadrupolar path combines field-gradient data with isotope spin and quadrupole moment. Testing just the element symbols or a rendered ellipsoid will not detect a wrong frequency scale.
