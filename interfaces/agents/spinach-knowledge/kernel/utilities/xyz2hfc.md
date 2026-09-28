# kernel/utilities/xyz2hfc.m

- Signature: `A=xyz2hfc(exyz,nxyz,isotope)`

## Purpose

Converts point electron and nuclear coordinates into a dipolar hyperfine interaction tensor.

## Parameters / inputs

- `exyz`: Electron Cartesian coordinates, a real `1x3` row vector in Angstrom.
- `nxyz`: Nuclear Cartesian coordinates, a real `1x3` row vector in Angstrom; it must differ from `exyz`.
- `isotope`: Isotope specification as a character string, for example `'13C'`.

## Numerical / algorithmic content

The displacement from electron to nucleus is `r=nxyz-exyz`. The function obtains the nuclear magnetogyric ratio as `gamma_n=spin(isotope)` and computes the dipolar matrix `D=3*(r'*r)/norm(r,2)^5-eye(3)/norm(r,2)^3`. It returns `A=C*D`, where `C=10^4*gamma_n*hbar*mu0/(4*pi*(1e-10)^3)`, `hbar=1.054571730e-34`, and `mu0=4*pi*1e-7`.

## Outputs

- `A`: Hyperfine coupling tensor in Gauss. Gauss units are used because the coupling does not depend on the electron `g`-tensor. The returned tensor enters the spin Hamiltonian as `S*A*I`; it is not scaled by the number of unpaired electrons because the electron spin operator already carries that magnitude.

https://spindynamics.org/wiki/index.php?title=xyz2hfc.m