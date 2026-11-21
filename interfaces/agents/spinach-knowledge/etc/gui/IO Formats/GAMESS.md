# etc/gui/IO Formats/GAMESS.cpp

[Implementation](https://github.com/IlyaKuprov/Spinach/blob/main/etc/gui/IO%20Formats/GAMESS.cpp)

The equilibrium-geometry marker identifies the Cartesian table, already in angstroms. Before an equilibrium frame is available the reader can use an initial ATOM/ATOMIC/COORDINATES (BOHR) table, multiplying by 0.529177249. An equilibrium geometry later in the file updates the atoms. This is a geometry import, not a general GAMESS property reader. A successful import does not imply that coupling or shielding tensors were read.
