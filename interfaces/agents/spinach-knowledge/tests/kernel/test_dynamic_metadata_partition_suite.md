# tests/kernel/test_dynamic_metadata_partition_suite.m

**Source:** [https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_dynamic_metadata_partition_suite.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_dynamic_metadata_partition_suite.m)

## Purpose

Regression test for deterministic metadata, hashing, and partition helper functions in Spinach. The suite verifies hashing stability, duplicate-row removal, parallel-state metadata, transfer matrix recovery, graph components, and safe partition exits.

## Behaviour

The test announces its target with `TESTING: Metadata and partition utilities` and initialises a regression test result via `new_test_result` for `kernel/dynamic_metadata_partition_suite`, described as `Metadata, hashing, and partition utilities`, with the requirement that `small metadata and partition helpers must preserve stable identity and exact graph or algebraic behaviour.`

The suite then performs the following checks:

- **Parallel-state metadata (client):** calls `poolsize()` and requires the result to be numeric, scalar, non-negative, and integer-valued (`mod(pool_count,1)==0`); calls `isworkernode()` and requires it to be false on the client, since the validation driver runs on the MATLAB client rather than a parallel worker.
- **MD5 hash stability and sensitivity:** computes `md5_hash({[1 2 3],'abc'})` twice and requires identical 32-character hexadecimal strings (`isstrprop(hash_a,'xdigit')`); requires `md5_hash({[1 2 4],'abc'})` to differ, so changing the serialised object contents changes the hash; requires `md5_hash(eye(2))` and `md5_hash(speye(2))` to differ, since full and sparse matrices are distinct MATLAB objects.
- **Duplicate-row removal:** builds the sparse matrix `A=sparse([1 0 2;1 0 2;0 3 0;1 0 2;4 0 0])` and checks `unihash(A)` against `sparse([1 0 2;0 3 0;4 0 0])` with tolerances `1e-15` (absolute and relative), verifying that the first occurrence of each unique sparse row is kept in stable order.
- **Transfer matrix recovery:** with `T_ref=[2 1;0 -1]`, `amp_inps=[1 0 1;0 1 1]`, and `amp_outs=T_ref*amp_inps`, checks `transfermat(amp_inps,amp_outs)` against `T_ref` with tolerances `1e-14` (absolute and relative), verifying that linearly complete input-output samples recover the exact linear transfer matrix.
- **Strongly connected components:** with `G=logical([1 1 0;1 1 0;0 0 1])`, checks `scomponents(G)` so that `sci(1)==sci(2)`, `sci(3)~=sci(1)`, and `numel(unique(sci))==2`; nodes one and two are mutually reachable and node three is a separate component.
- **Path tracing disabled exit:** sets `spin_system.sys.output='hush'` and `spin_system.sys.disable={'pt'}`, then calls `path_trace(spin_system,speye(3),[])`; the result must be a scalar cell whose first element equals `1`, i.e. a unit projector placeholder returned without graph partition work.
- **Zero-track elimination disabled exit:** sets `spin_system.bas.formalism='sphten-liouv'` and `spin_system.sys.enable={}`, then calls `zte(spin_system,speye(3),[1;0;0])`; the result must equal `1`, i.e. a unit projector placeholder returned without Krylov propagation.

## Inputs and outputs

Syntax:

```matlab
result=test_dynamic_metadata_partition_suite()
```

The function takes no inputs. It returns `result`, a regression test result structure with explanatory messages accumulated by `test_true` and `test_close` checks.

## References

- [Spinach source: tests/kernel/test_dynamic_metadata_partition_suite.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_dynamic_metadata_partition_suite.m)

The suite also verifies that explicit ZTE enablement removes three empty tracks from a four-coordinate invariant system, and that paranoia overrides the enablement.
