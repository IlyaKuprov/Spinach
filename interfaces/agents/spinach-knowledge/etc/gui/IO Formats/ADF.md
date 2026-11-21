# etc/gui/IO Formats/ADF.cpp

[Implementation](https://github.com/IlyaKuprov/Spinach/blob/main/etc/gui/IO%20Formats/ADF.cpp)

ADF geometry, shielding and scalar couplings are recognized by distinct output headings. The shielding tensor belongs to one atom and is stored in ppm; the total scalar J belongs to a pair and is stored in Hz. A pair printed in both directions should not double the number of couplings. For repeated geometry/property sections check which atoms the retained tensor references, not only whether the file opened.
