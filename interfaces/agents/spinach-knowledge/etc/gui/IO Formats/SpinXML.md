# etc/gui/IO Formats/SpinXML.cpp

[Implementation](https://github.com/IlyaKuprov/Spinach/blob/main/etc/gui/IO%20Formats/SpinXML.cpp)

Use SpinXML for the editable spin-system document. The reader maps external spin and reference-frame IDs to model objects; references need not be dense, one-based or listed in parent-before-child order. Missing spin_b means a one-spin interaction, not a reference to the first nucleus. An omitted interaction frame defaults to the lab frame.

Tensor input may be given as a matrix or in principal-value representations with an orientation. On round trip compare actual atom pairs, lab-frame matrices, units and frame relationships; textual identity or unchanged external numbering is not the criterion. Molecular susceptibility is separate from the interaction collection and is not serialized by this writer. Preserve the Gaussian source or a susceptibility-bearing export if that property is needed later.
