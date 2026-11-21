# etc/gui/IO Formats/Simpson.cpp

[Implementation](https://github.com/IlyaKuprov/Spinach/blob/main/etc/gui/IO%20Formats/Simpson.cpp)

The SIMPSON writer builds a spinsys block with channels and nuclei. Channel names identify isotope types, whereas interaction indices refer to individual spins in the exported order. When choosing channels or sorting atoms, compare the emitted pair indices against that order. Chemical shifts, scalar couplings, quadrupolar terms and bond-derived dipolar terms have different syntax; a bond drawn in the GUI is not itself a scalar J. Inspect the final text before using it as an experiment input.
