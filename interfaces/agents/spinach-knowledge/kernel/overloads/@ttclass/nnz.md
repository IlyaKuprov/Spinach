# kernel/overloads/@ttclass/nnz.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/overloads/@ttclass/nnz.m`
- Signature: `answer=nnz(ttrain)`
- Total lines: 33

## Purpose

Counts non-zero elements in all cores of a tensor train. Syntax: answer=nnz(ttrain)

## Physical / mathematical content

- Tensor-train linear algebra. These files implement compressed high-dimensional operators and AMEn/SVD-based algebra in tensor-train format.

## Numerical / algorithmic content


## Code-derived implementation details

### Comment-guided execution stages

- Lines 20-21: Count the non-zeros; implemented by `answer=sum(sum(cellfun(@nnz,ttrain.cores)))`.

### Key state/data transformations

- Lines 21: computes `answer` using `answer=sum(sum(cellfun(@nnz,ttrain.cores)))`.

## Parameters / inputs

- ttrain -tensor train object

## Outputs

- answer -number of non-zero elements in all tensor train cores

## Implementation structure

- Counts non-zero elements in all cores of a tensor train. Syntax:
- answer=nnz(ttrain)
- ttrain -tensor train object
- answer -number of non-zero elements in all tensor train cores
- Count the non-zeros
- Men have always been and forever would remain silly victims of lies and
- self-deceit in politics, until they learn to see, behind any moral, re-
- ligious, political or social statements, proclamations and promises the
- interests of specific social classes.
- Vladimir Lenin
- #NGRUM

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `cellfun()`.
