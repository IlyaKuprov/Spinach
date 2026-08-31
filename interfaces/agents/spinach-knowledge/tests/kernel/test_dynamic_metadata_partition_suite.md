# tests/kernel/test_dynamic_metadata_partition_suite.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_dynamic_metadata_partition_suite.m`
- Signature: `result=test_dynamic_metadata_partition_suite()`
- Total lines: 78

## Purpose

Tests deterministic metadata, hashing, and partition helpers. Syntax: result=test_dynamic_metadata_partition_suite()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.
- Propagation is accelerated with a Krylov-subspace method, replacing direct matrix exponentiation by projection into a much smaller Arnoldi/Lanczos-type subspace.

## Numerical / algorithmic content

- A Krylov-subspace or Arnoldi construction is used to avoid forming or exponentiating very large dense propagators directly.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 16-17: Announce the test target; implemented by `fprintf('TESTING: Metadata and partition utilities\n')`.
- Lines 19-22: State the metadata and partition target of the test; implemented by `result=new_test_result('kernel/dynamic_metadata_partition_suite', 'Metadata, hashing, and partition utilities', 'small metadata and partition helpers must preserve stabl…`.
- Lines 24-25: Check parallel-state metadata on the MATLAB client; implemented by `pool_count=poolsize()`.
- Lines 33-34: Check MD5 hash stability and object-type sensitivity; implemented by `hash_a=md5_hash({[1 2 3],'abc'})`.
- Lines 44-45: Check stable duplicate-row removal through hash-table identity; implemented by `A=sparse([1 0 2;1 0 2;0 3 0;1 0 2;4 0 0])`.
- Lines 49-50: Check least-squares transfer matrix recovery from overdetermined samples; implemented by `T_ref=[2 1;0 -1]`.
- Lines 56-57: Check strongly connected components on a two-component directed graph; implemented by `G=logical([1 1 0;1 1 0;0 0 1])`.
- Lines 62-63: Check path tracing disabled exit without graph partition work; implemented by `spin_system.sys.output='hush'`.
- Lines 69-70: Check zero-track elimination disabled exit without Krylov propagation; implemented by `spin_system.bas.formalism='sphten-liouv'`.

### Key state/data transformations

- Lines 20-22: computes `result` using `result=new_test_result('kernel/dynamic_metadata_partition_suite', 'Metadata, hashing, and partition utilities', 'small metadata and partition helpers must preserve stabl…`.
- Lines 25: computes `pool_count` using `pool_count=poolsize()`.
- Lines 34: computes `hash_a` using `hash_a=md5_hash({[1 2 3],'abc'})`.
- Lines 35: computes `hash_b` using `hash_b=md5_hash({[1 2 3],'abc'})`.
- Lines 36: computes `hash_c` using `hash_c=md5_hash({[1 2 4],'abc'})`.
- Lines 45: computes `A` using `A=sparse([1 0 2;1 0 2;0 3 0;1 0 2;4 0 0])`.
- Lines 50: computes `T_ref` using `T_ref=[2 1;0 -1]`.
- Lines 51: computes `amp_inps` using `amp_inps=[1 0 1;0 1 1]`.
- Lines 52: computes `amp_outs` using `amp_outs=T_ref*amp_inps`.
- Lines 57: computes `G` using `G=logical([1 1 0;1 1 0;0 0 1])`.
- Lines 58: computes `sci` using `sci=scomponents(G)`.
- Lines 63: computes `spin_system.sys.output` using `spin_system.sys.output='hush'`.
- Lines 64: computes `spin_system.sys.disable` using `spin_system.sys.disable={'pt'}`.
- Lines 65: computes `projectors` using `projectors=path_trace(spin_system,speye(3),[])`.
- Lines 70: computes `spin_system.bas.formalism` using `spin_system.bas.formalism='sphten-liouv'`.
- Lines 72: computes `projector` using `projector=zte(spin_system,speye(3),[1;0;0])`.

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

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `poolsize()`, `test_true()`, `isscalar()`, `isworkernode()`, `md5_hash()`, `strcmp()`, `all()`, `isstrprop()`, `speye()`, `test_close()`, `unihash()`, `transfermat()`, `logical()`, `scomponents()`, `sci()`.
