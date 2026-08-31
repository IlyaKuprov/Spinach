# etc/textbook/lorentz.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/etc/textbook/lorentz.m`
- Signature: `[J,K,Kil]=lorentz(L)`
- Total lines: 83

## Purpose

The (L,0)(+)(0,L) irreducible matrix representation of the Lorentz group with inversion. Syntax: [J,K,Kil]=lorentz(L)

## Physical / mathematical content

- This file belongs to the `etc` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 25-26: Check consistency; implemented by `grumble(L)`.
- Lines 28-29: Dimension and Pauli blocks; implemented by `D=2*L+1; s=pauli(D)`.
- Lines 31-32: Rotation and boost generators; implemented by `J.x=full([s.x zeros(D,D); zeros(D,D) s.x])`.
- Lines 39-40: Collect generators; implemented by `gens={J.x J.y J.z K.x K.y K.z}`.
- Lines 42-43: Killing form is expensive; implemented by `if nargout>2`.
- Lines 45-46: Preallocate Killing form; implemented by `Kil=zeros(numel(gens),numel(gens),'like',1i)`.
- Lines 48-49: Loop over generators; implemented by `for n=1:numel(gens)`.
- Lines 51-52: Get the first operator adjoint; implemented by `AdA=kron(eye(2*D,2*D),gens{n})-kron(transpose(gens{n}),eye(2*D,2*D))`.
- Lines 54-55: Loop over generators; implemented by `for k=1:numel(gens)`.
- Lines 57-58: Get the second operator adjoint; implemented by `AdB=kron(eye(2*D,2*D),gens{k})-kron(transpose(gens{k}),eye(2*D,2*D))`.
- Lines 60-61: Get the Killing form element; implemented by `Kil(n,k)=trace(AdA*AdB)`.

### Control flow inferred from the code

- Line 43: conditional branch on `nargout>2`.
- Line 49: `for` loop over `n=1:numel(gens)`.
- Line 55: `for` loop over `k=1:numel(gens)`.

### Key state/data transformations

- Lines 29: computes `D` using `D=2*L+1; s=pauli(D)`.
- Lines 32: computes `J.x` using `J.x=full([s.x zeros(D,D); zeros(D,D) s.x])`.
- Lines 33: computes `J.y` using `J.y=full([s.y zeros(D,D); zeros(D,D) s.y])`.
- Lines 34: computes `J.z` using `J.z=full([s.z zeros(D,D); zeros(D,D) s.z])`.
- Lines 35: computes `K.x` using `K.x=full([1i*s.x zeros(D,D); zeros(D,D) -1i*s.x])`.
- Lines 36: computes `K.y` using `K.y=full([1i*s.y zeros(D,D); zeros(D,D) -1i*s.y])`.
- Lines 37: computes `K.z` using `K.z=full([1i*s.z zeros(D,D); zeros(D,D) -1i*s.z])`.
- Lines 40: computes `gens` using `gens={J.x J.y J.z K.x K.y K.z}`.
- Lines 46: computes `Kil` using `Kil=zeros(numel(gens),numel(gens),'like',1i)`.
- Lines 52: computes `AdA` using `AdA=kron(eye(2*D,2*D),gens{n})-kron(transpose(gens{n}),eye(2*D,2*D))`.
- Lines 58: computes `AdB` using `AdB=kron(eye(2*D,2*D),gens{k})-kron(transpose(gens{k}),eye(2*D,2*D))`.
- Lines 61: computes `Kil(n,k)` using `Kil(n,k)=trace(AdA*AdB)`.

### Local helper functions

- Line 72: `grumble()` — `function grumble(L)`. The FBI is the only organization on Earth complaining that computer security is too good.
  - Representative operation: `if (~isnumeric(L))||(~isscalar(L))||(~isreal(L))|| (mod(2*L,1)~=0)||(L<0.5)`.
  - Representative operation: `(mod(2*L,1)~=0)||(L<0.5)`.

## Parameters / inputs

- L -irreducible representation
- rank, e.g. 1/2

## Outputs

- J -three rotation generators
- K -three boost generators
- Kil -Killing form (expensive)

## Implementation structure

- The (L,0)(+)(0,L) irreducible matrix representation of the
- Lorentz group with inversion. Syntax:
- [J,K,Kil]=lorentz(L)
- L - irreducible representation
- rank, e.g. 1/2
- J - three rotation generators
- K - three boost generators
- Kil - Killing form (expensive)
- Check consistency
- Dimension and Pauli blocks
- Rotation and boost generators
- Collect generators

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `pauli()`, `transpose()`, `Kil()`, `isscalar()`.
