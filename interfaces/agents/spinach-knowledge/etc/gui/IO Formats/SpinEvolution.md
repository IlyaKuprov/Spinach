# etc/gui/IO Formats/SpinEvolution.cpp

[Implementation](https://github.com/IlyaKuprov/Spinach/blob/main/etc/gui/IO%20Formats/SpinEvolution.cpp)

The writer generates a SpinEvolution text input with channel labels, a nuclei list and interaction and pulse sections. Filtering channels changes the set of exported atoms, so validate shift, scalar-coupling and quadrupolar indices against the output rather than the grid row numbers. The template includes spectrometer(MHz) 500; an example field or pulse sequence in generated text is not information measured in the imported structure.
