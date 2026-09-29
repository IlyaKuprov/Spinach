# examples/fundamentals/operator_tests/commutation_7.m

- Signature: `commutation_7()`

## Purpose

Checks that finite-dimensional operators can be reconstructed consistently from several Spinach operator bases and expansion routines. The examples cover spin Zeeman level projectors, bosonic number-level projectors and ladder-operator products, an arbitrary matrix in a single-transition basis, and an arbitrary bosonic matrix in a monomial basis.

## Mathematical content

The question is whether the coefficients returned by the expansion helpers reproduce the operator represented by the selected basis elements. This is a basis-conversion check, not a simulation of a measured spin or boson response. In particular, the bosonic product cases construct explicit products of the finite-matrix creation and annihilation operators and compare them with an irreducible spherical tensor (IST) expansion.

## Callable context and model

Call the zero-input MATLAB function `commutation_7()` from a Spinach checkout with the project functions on the MATLAB path. It has no input or returned output; it prints a success message for each group and raises an error when a stated comparison exceeds its threshold. The test uses a spin multiplicity of 5 and a bosonic truncation of 6 levels. Its matrix examples are complex random matrices of dimensions 7 and 6, respectively.

## Checks encoded in the source

- For each of the 5 Zeeman levels, forms a diagonal projector, obtains coefficients with `enlev2ist(5,lvl_num,'S')`, and reconstructs it from `irr_sph_ten(5)`. The Frobenius-norm difference must be at most `1e-10`.
- For each of the 6 bosonic levels, forms a diagonal projector, expands it with `enlev2bm(6,lvl_num)`, and reconstructs it from `boson_mono(6)`. The Frobenius-norm difference must be at most `1e-10`.
- Builds the three literal ladder strings `CA`, `ACCA`, and `CCAAA` using `W.c` for each `C` and `W.a` for each `A`, where `W=weyl(6)`. It compares each product with its `bos2ist` expansion in `irr_sph_ten(6)`; the absolute Frobenius-norm difference must be at most `1e-10`.
- Expands a complex random 7-by-7 matrix in `sin_tran(7)`, using `hdot(E{n},A)` as each coefficient. The absolute Frobenius-norm difference must be at most `1e-10`.
- Expands a complex random 6-by-6 matrix using `oper2bm` and reconstructs it from `boson_mono(6)`. Its relative Frobenius-norm difference, divided by the original matrix norm, must be at most `1e-8`.

## Assumptions and limits

The finite dimensions and operator products are fixed by this function; the two matrix examples use unseeded random values. Thus the source defines concrete numerical tolerances but does not contain a fixed set of output values. Passing these checks would establish only the tested finite-dimensional reconstruction identities and conventions; it would not establish completeness for other dimensions or validate a physical model. The success strings in the source are conditional messages, not a result claimed here.

## Source

[`examples/fundamentals/operator_tests/commutation_7.m`](https://github.com/IlyaKuprov/Spinach/blob/main/examples/fundamentals/operator_tests/commutation_7.m)
