# kernel/overloads/@polyadic/ctranspose.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/overloads/@polyadic/ctranspose.m`
- Signature: `p=ctranspose(p)`
- Total lines: 42

## Purpose

Computes the Hermitian conjugate of a matrix in a polyadic representation. Syntax: p=ctranspose(p)

## Physical / mathematical content

- Polyadic tensor-product linear algebra. The emphasis is compressed operator representation, deferred algebra, and efficient Kronecker-structured manipulations.

## Numerical / algorithmic content


## Code-derived implementation details

### Comment-guided execution stages

- Lines 12-13: Conjugate-transpose every core; implemented by `for n=1:numel(p.cores)`.
- Lines 19-20: Process prefix and suffix; implemented by `new_prefix=fliplr(p.suffix)`.

### Control flow inferred from the code

- Line 13: `for` loop over `n=1:numel(p.cores)`.
- Line 14: `for` loop over `k=1:numel(p.cores{n})`.
- Line 21: `for` loop over `n=1:numel(new_prefix)`.
- Line 25: `for` loop over `n=1:numel(new_suffix)`.

### Key state/data transformations

- Lines 15: computes `p.cores{n}{k}` using `p.cores{n}{k}=ctranspose(p.cores{n}{k})`.
- Lines 20: computes `new_prefix` using `new_prefix=fliplr(p.suffix)`.
- Lines 22: computes `new_prefix{n}` using `new_prefix{n}=ctranspose(new_prefix{n})`.
- Lines 24: computes `new_suffix` using `new_suffix=fliplr(p.prefix)`.
- Lines 26: computes `new_suffix{n}` using `new_suffix{n}=ctranspose(new_suffix{n})`.
- Lines 28: computes `p.prefix` using `p.prefix=new_prefix; p.suffix=new_suffix`.

## Implementation structure

- Computes the Hermitian conjugate of a matrix in a polyadic
- representation. Syntax:
- p=ctranspose(p)
- Conjugate-transpose every core
- Process prefix and suffix
- Guys, stop wasting your time. There can be no such thing as
- a personal computer. Personal car, personal pension, perso-
- nal dacha -perhaps. Do you even know what a computer is? A
- computer is 100 square metres of building space, 25 service
- staff members, and 30 litres of ethanol a month!
- Nikolai Gorshkov, USSR Deputy Minister
- for Radioelectronics Industry, 1980

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `fliplr()`.
