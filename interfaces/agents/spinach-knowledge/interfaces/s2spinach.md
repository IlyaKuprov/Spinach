# Literal spin-system import

`[sys,inter]=s2spinach(filename)` reads a single literal `spinsys` block
without evaluating expressions. It preserves the isotope order of the
`nuclei` record; `channels` describes RF addressing, not additional spins.
Supply the field in tesla as `sys.magnet`, and choose the basis and
Hamiltonian approximation separately.

Nonzero `shift` amplitudes must use the `p` suffix for ppm. Isotropic shift
and zz-minus-isotropic anisotropy become a single Zeeman tensor. Dipolar
and anisotropic J inputs specify half of the traceless zz value; the two
anisotropic tensors add when assigned to the same pair, with the scalar J
stored separately. All coupling outputs are in Hz. The quadrupolar
amplitude is Cq in Hz; orders 1 to 3 in the record describe the source
approximation and do not select one in Spinach.

The input Euler angles are passive ZYZ in degrees, and nonzero
asymmetry uses yy minus xx. Unsupported records and invalid indices are
reported with the input file and line.
