# tests/kernel/test_vstore.m

Tests the current parallel pool ValueStore save/merge contract using a parallel pool obtained through Spinach; the fixture requests one worker when creating a new pool. The snapshot includes a complex dense matrix, a sparse identity matrix, a nested structure, and a string. Checks cover saving without mutation, exact value round-trips, sparse representation, preservation of current-only keys, replacement of matching values, the union of key sets, repeated loading, an empty snapshot, and rejection of a malformed snapshot without changing the store.
