# etc/gui/IO Formats/XYZ.cpp

[Implementation](https://github.com/IlyaKuprov/Spinach/blob/main/etc/gui/IO%20Formats/XYZ.cpp)

Use XYZ for geometry, not spin tensors. Each frame starts with a non-negative integer atom count and a comment line, followed by exactly that many four-field element/x/y/z rows. Even a comment that looks like an atom is skipped. Coordinates have scale one in the model's angstrom convention. Repeated frames update the geometry; the final refresh calculates display bonds. There are no isotope-specific coupling or shielding sections in this reader.
