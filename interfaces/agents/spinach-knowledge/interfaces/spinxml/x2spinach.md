# interfaces/spinxml/x2spinach.m

Source: [interfaces/spinxml/x2spinach.m](https://github.com/IlyaKuprov/Spinach/blob/main/interfaces/spinxml/x2spinach.m)

Signature: `[sys,inter] = x2spinach(filename,shielding_refs)`

Reads a SpinXML file and constructs Spinach spin-system and interaction input fields.

## Input contract

- filename must be a non-empty character string naming an existing file.
- shielding_refs is a cell array of two-item cells: an isotope string and its absolute shielding reference. Isotope entries must be unique. For example, {{'1H',31.5},{'13C',189.7}}. If the XML supplies chemical shifts rather than absolute shieldings, use an empty cell array.
- The source warns that the file must already have passed validation against the [SpinXML schema](http://spindynamics.org/SpinXML.php). This routine calls parsexml(filename); it does not itself run schema validation.

## Spin and interaction data

The first pass reads each spin element's ID and isotope, and its optional label and Cartesian coordinate attributes x, y, and z. Spins are sorted by ID; the IDs must then be exactly the sequence 1 through the number of spins. The coordinate values are copied into the output without a unit conversion in this function.

A second pass reads interaction elements, including their kind, spin IDs, units, labels, references, tensor/scalar specification, and optional orientation. The value forms handled in the source are tensor, scalar, eigenvalues, aniso_asymm, axiality_rhombicity, and span_skew. Recognised interaction kinds are hfc, shielding, shift, dipolar, quadrupolar, jcoupling, gtensor, zfs, exchange, and spinrotation.

For absolute shielding, the supplied isotope reference is used to form the chemical-shift tensor as reference * I - shielding; the interaction must be nuclear and the reference for that isotope must be supplied. Chemical shifts and shieldings require ppm; gtensor requires bohr. The frequency-valued hfc, dipolar, quadrupolar, jcoupling, zfs, exchange, and spinrotation branches accept Hz, kHz, MHz, or GHz; hfc also accepts gauss. The source rejects unsupported units and incompatible spin assignments.

## Returned fields and dimensions

For N spins, the routine populates sys.isotopes and sys.labels as one-entry-per-spin cell arrays, and inter.coordinates as an N×1 cell array containing each coordinate vector or an empty value. Interaction storage is cell-based: inter.zeeman.matrix, inter.zeeman.reference, and inter.zeeman.label are 1×N; inter.coupling.matrix and inter.coupling.label are N×N; and inter.spinrot.matrix and inter.spinrot.label are 1×N. Tensor-valued entries are stored as 3×3 matrices. These are interaction input tensors, not Hilbert-space operators or state vectors.

The source also retains the SpinXML warning about a possible double count when both spins in a dipolar interaction have coordinates. It does not define a cache or remote-file retrieval layer.

## Reference

[Spin Dynamics Wiki: x2spinach.m](https://spindynamics.org/wiki/index.php?title=x2spinach.m)
