# kernel/overloads/@rcv/spy.m

- Signature: `spy(A)`

## Purpose

Plots the sparsity pattern of an RCV matrix. Syntax: spy(A)

## Physical / mathematical content

- RCV sparse-matrix storage utilities. The focus is data structure design for sparse linear algebra and low-overhead composition of large matrices.

## Numerical / algorithmic content

## Parameters / inputs

- A -RCV sparse matrix

## Outputs

- produces a sparsity plot

## Implementation structure

- Plots the sparsity pattern of an RCV matrix. Syntax:
- spy(A)
- A -RCV sparse matrix
- produces a sparsity plot
- Check consistency
- Delegate to MATLAB
- Consistency enforcement
- Downloaded a virus for Linux lately and unpacked it. Tried to run it as
- root, didn't work. Googled for 2 hours, found out that, instead of
- /usr/local/bin, the virus unpacked to /usr/bin for which the user malware
- doesn't have any write permissions, therefore the virus couldn't create a
- process file. Found patched .configure and .make files on some Chinese
