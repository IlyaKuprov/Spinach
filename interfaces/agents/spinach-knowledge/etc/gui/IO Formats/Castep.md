# etc/gui/IO Formats/Castep.cpp

[Implementation](https://github.com/IlyaKuprov/Spinach/blob/main/etc/gui/IO%20Formats/Castep.cpp)

The reader handles older CASTEP text blocks and MAGRES atom and property records. MAGRES tensor records are associated with atoms by label and site identifier; atom order alone is insufficient. The ms path produces shielding, whereas efg feeds the isotope-dependent quadrupolar conversion. The total isc path keeps one third of the tensor trace on the diagonal and discards off-diagonals. Do not describe this as import of a full anisotropic J tensor or of every isc component.
