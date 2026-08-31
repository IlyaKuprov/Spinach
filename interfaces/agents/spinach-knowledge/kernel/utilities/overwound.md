# kernel/utilities/overwound.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/overwound.m`
- Signature: `overwound(rho,spc_dim,spn_dim)`
- Total lines: 140

## Purpose

Checks if a Fokker-Planck state vector has any spatial frequencies that its spatial grid is dangerously close to misrepresenting due to insufficient point count. Syntax: overwound(rho,spc_dim,spn_dim)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.
- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 33-34: Check consistency; implemented by `grumble(rho,spc_dim,spn_dim)`.
- Lines 36-37: Find out the stack size; implemented by `stack_size=size(rho,2)`.
- Lines 39-40: Fold the state vector; implemented by `rho=reshape(rho,[spn_dim spc_dim stack_size])`.
- Lines 42-43: Get spin dimensions out of the way; implemented by `rho=permute(rho,[2 3 4 1 5])`.
- Lines 45-46: Check X dimension; implemented by `if spc_dim(1)>1`.
- Lines 48-49: Few points are not good; implemented by `if spc_dim(1)<10`.
- Lines 53-54: Run Fourier transform; implemented by `rho_kx=fftshift(fft(rho,[],1),1)`.
- Lines 56-57: Pull X dimension forward; implemented by `rho_kx=permute(rho_kx,[1 2 3 4 5])`.
- Lines 61-62: Draw the plot; implemented by `freq_axis=linspace(-1,1,size(rho_kx,1))`.
- Lines 69-70: Check Y dimension; implemented by `if spc_dim(2)>1`.
- Lines 72-73: Few points are not good; implemented by `if spc_dim(2)<10`.
- Lines 77-78: Run Fourier transform; implemented by `rho_ky=fftshift(fft(rho,[],2),2)`.
- Lines 80-81: Pull Y dimension forward; implemented by `rho_ky=permute(rho_ky,[2 1 3 4 5])`.
- Lines 85-86: Draw the plot; implemented by `freq_axis=linspace(-1,1,size(rho_ky,1))`.
- Lines 93-94: Check X dimension; implemented by `if spc_dim(3)>1`.
- Lines 96-97: Few points are not good; implemented by `if spc_dim(3)<10`.
- Lines 101-102: Run Fourier transform; implemented by `rho_kz=fftshift(fft(rho,[],3),3)`.
- Lines 104-105: Pull X dimension forward; implemented by `rho_kz=permute(rho_kz,[3 2 1 4 5])`.

### Control flow inferred from the code

- Line 46: conditional branch on `spc_dim(1)>1`.
- Line 49: conditional branch on `spc_dim(1)<10`.
- Line 70: conditional branch on `spc_dim(2)>1`.
- Line 73: conditional branch on `spc_dim(2)<10`.
- Line 94: conditional branch on `spc_dim(3)>1`.
- Line 97: conditional branch on `spc_dim(3)<10`.

### Key state/data transformations

- Lines 37: computes `stack_size` using `stack_size=size(rho,2)`.
- Lines 40: computes `rho` using `rho=reshape(rho,[spn_dim spc_dim stack_size])`.
- Lines 54: computes `rho_kx` using `rho_kx=fftshift(fft(rho,[],1),1)`.
- Lines 62: computes `freq_axis` using `freq_axis=linspace(-1,1,size(rho_kx,1))`.
- Lines 78: computes `rho_ky` using `rho_ky=fftshift(fft(rho,[],2),2)`.
- Lines 102: computes `rho_kz` using `rho_kz=fftshift(fft(rho,[],3),3)`.

### Local helper functions

- Line 120: `grumble()` — `function grumble(rho,spc_dim,spn_dim)`.
  - Representative operation: `if (~isnumeric(rho))`.
  - Representative operation: `error('rho must be a numeric array.')`.

## Parameters / inputs

- rho -Fokker-Planck state vector or
- a bookshelf stack thereof
- spc_dim -spatial dimensions of the Fokker-
- Planck problem, [X Y Z]
- spn_dim -spin dimension of the Fokker-
- Planck problem
- Output:
- figures and diagnostic messages to the console
- Note: if you are running spatial dynamics, such as diffusion and
- and flow, with finite difference derivative operators, the
- spatial grid point count shuld be set to several times the
- minimum Nyquist value.

## Implementation structure

- Checks if a Fokker-Planck state vector has any spatial frequencies
- that its spatial grid is dangerously close to misrepresenting due
- to insufficient point count. Syntax:
- overwound(rho,spc_dim,spn_dim)
- rho -Fokker-Planck state vector or
- a bookshelf stack thereof
- spc_dim -spatial dimensions of the Fokker-
- Planck problem, [X Y Z]
- spn_dim -spin dimension of the Fokker-
- Planck problem
- Output:
- figures and diagnostic messages to the console

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `spc_dim()`, `fftshift()`, `kfigure()`, `kxlabel()`, `kylabel()`, `isvector()`, `any()`, `isscalar()`.
