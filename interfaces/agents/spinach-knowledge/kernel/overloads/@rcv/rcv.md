# kernel/overloads/@rcv/rcv.m

- Signature: `obj=rcv(varargin)`

## Purpose

Constructs an RCV sparse-matrix object from a matrix, dimensions, or explicit row, column, and value arrays.

## Mathematical content

The object stores matrix entries as parallel row-index, column-index, and value arrays, together with the matrix dimensions. For a matrix input, the constructor obtains these arrays from its nonzero entries.

## Numerical / algorithmic content

The one-input form returns an existing RCV object unchanged or converts a matrix by extracting its nonzeros. The two-dimension form creates an empty CPU-resident object. The five-input form casts indices and dimensions to `int64`, values to `double`, and marks the object GPU-resident if any entry-array input is a `gpuArray`; in that case the stored arrays are uploaded to the GPU. It validates argument count, dimensions, array lengths, and index bounds.

## Parameters / inputs

- M -a Matlab matrix
- dim1 -number of rows
- dim2 -number of columns
- R -row indices of non-zero entries
- C -column indices of non-zero entries
- V -values corresponding to entries in R and C

## Outputs

- obj -an RCV sparse matrix object

## Implementation structure

- Validate inputs for the selected one-, two-, or five-argument form.
- For a matrix input, preserve its GPU location, record its dimensions, extract nonzero entries, and store row and column indices as `int64` and values as `double`.
- For explicit arrays, store the supplied entries and dimensions using those types; upload the entry arrays when any of them is GPU-resident.
- The class also reports `true` for `isnumeric`, `ismatrix`, and `isfloat`.
