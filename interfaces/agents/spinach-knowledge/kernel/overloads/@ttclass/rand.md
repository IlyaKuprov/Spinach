# kernel/overloads/@ttclass/rand.m

- Signature: `tt=rand(tt,ttrank)`

## Purpose

Generates a tensor train with random cores, retaining the supplied train's physical dimensions and using the requested bond rank for internal bonds (except that a one-core train has boundary ranks 1).

## Physical / mathematical content

The generated cores have the same physical index dimensions as the input. The result has coefficient 1 and tolerance 0.

## Numerical / algorithmic content

Each core is filled with MATLAB `rand`. For multiple cores, the first and last ranks are 1 and each internal bond has rank `ttrank`.

## Parameters / inputs

- tt - a tensor train object
- ttrank - bond rank, a positive real integer

## Outputs

- tt - a tensor train object

## Implementation structure

- Validate the object and rank input.
- Read the physical sizes, allocate the cores, and fill them with random values.
- Set the coefficient to 1 and tolerance to 0.
