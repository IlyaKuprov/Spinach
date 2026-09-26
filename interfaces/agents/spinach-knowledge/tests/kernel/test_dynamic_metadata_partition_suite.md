# tests/kernel/test_dynamic_metadata_partition_suite.m

- Signature: `result=test_dynamic_metadata_partition_suite()`

## Purpose

Tests deterministic metadata, hashing, and partition helpers. Syntax: result=test_dynamic_metadata_partition_suite()

## Physical / mathematical content

- Propagation is accelerated with a Krylov-subspace method, replacing direct matrix exponentiation by projection into a much smaller Arnoldi/Lanczos-type subspace.

## Numerical / algorithmic content

- A Krylov-subspace or Arnoldi construction is used to avoid forming or exponentiating very large dense propagators directly.

## Outputs

- result -regression test result with explanatory messages
- The test checks hashing stability, duplicate-row removal, parallel-state
- metadata, transfer matrices, graph components, and safe partition exits.

## Implementation structure

- Tests deterministic metadata, hashing, and partition helpers. Syntax:
- result=test_dynamic_metadata_partition_suite()
- result -regression test result with explanatory messages
- The test checks hashing stability, duplicate-row removal, parallel-state
- metadata, transfer matrices, graph components, and safe partition exits.
- Announce the test target
- State the metadata and partition target of the test
- Check parallel-state metadata on the MATLAB client
- Check MD5 hash stability and object-type sensitivity
- Check stable duplicate-row removal through hash-table identity
- Check least-squares transfer matrix recovery from overdetermined samples
- Check strongly connected components on a two-component directed graph
