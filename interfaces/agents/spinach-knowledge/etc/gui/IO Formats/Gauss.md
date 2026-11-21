# etc/gui/IO Formats/Gauss.cpp

[Implementation](https://github.com/IlyaKuprov/Spinach/blob/main/etc/gui/IO%20Formats/Gauss.cpp)

The Gaussian reader recognizes separate geometry, isotope and magnetic-property sections. These need not all be present. Geometry comes from orientation tables; property blocks can supply shieldings, scalar J, g, hyperfine and quadrupolar interactions. For anisotropic hyperfine output the Baa/Bbb/Bcc rows are not interchangeable: the middle row includes an atom identifier and isotope. Combined contact and anisotropic data are used to form the electron–nucleus interaction.

Molecular magnetic susceptibility is read from the Gaussian cgs-ppm block into the model's separate susceptibility tensor. Do not invent an electron or nucleus to hold it, or infer failure from an empty interaction grid. Check its display separately and inspect the susceptibility fields in a Spinach export.
