# tests/kernel/test_perm_group_suite.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_perm_group_suite.m`
- Signature: `result=test_perm_group_suite()`
- Total lines: 142

## Purpose

Tests permutation group database metadata. Syntax: result=test_perm_group_suite()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.
- The relevant state manifold is the singlet/triplet decomposition, where permutation symmetry controls selection rules, relaxation susceptibility, and convertibility to ordinary magnetisation.

## Numerical / algorithmic content


## Code-derived implementation details

### Comment-guided execution stages

- Lines 17-18: Announce the test target; implemented by `fprintf('TESTING: Permutation group metadata\n')`.
- Lines 20-23: State the utility target of the test; implemented by `result=new_test_result('kernel/perm_group_suite', 'Permutation group metadata', 'real-valued-character Abelian permutation subgroup tables must be exact and internally c…`.
- Lines 25-26: Check the S6 Abelian subgroup table; implemented by `G=perm_group('S6A')`.
- Lines 42-43: Check the S8 Abelian subgroup table; implemented by `G=perm_group('S8A')`.
- Lines 69-70: Check structural consistency of all real-valued-character Abelian options; implemented by `group_names={'S4A','S6A','S8A'}`.
- Lines 75-76: Get the group under test; implemented by `G=perm_group(group_names{n})`.
- Lines 78-81: Check scalar group metadata; implemented by `result=test_true(result,[group_names{n} ' scalar metadata'], (G.order==group_orders(n))&&(G.nclasses==group_orders(n))&&(G.n_irreps==group_orders(n)), 'Abelian real-char…`.
- Lines 86-88: Check class and irrep dimensions; implemented by `result=test_close(result,[group_names{n} ' singleton classes'],G.class_sizes,ones(1,group_orders(n)),1e-15,1e-15, 'Abelian subgroup conjugacy classes are singletons')`.
- Lines 92-95: Check character reality and orthogonality; implemented by `result=test_true(result,[group_names{n} ' real characters'], isreal(G.class_characters)&&all(abs(G.class_characters(:))==1), 'real-valued-character Abelian subgroup char…`.
- Lines 100-101: Initialise multiplication checks; implemented by `group_closed=true; group_comm=true; group_invol=true; group_perm=true`.
- Lines 104-105: Check every permutation and every product; implemented by `for k=1:G.order`.
- Lines 117-119: Check permutation-table structure; implemented by `result=test_true(result,[group_names{n} ' permutation rows'],group_perm, 'all group elements must be valid permutation rows')`.
- Lines 127-128: Check maximality inside the ambient symmetric group; implemented by `central_count=0; all_perms=perms(identity)`.

### Control flow inferred from the code

- Line 73: `for` loop over `n=1:numel(group_names)`.
- Line 105: `for` loop over `k=1:G.order`.
- Line 109: `for` loop over `m=1:G.order`.
- Line 129: `for` loop over `k=1:size(all_perms,1)`.
- Line 131: `for` loop over `m=1:G.order`.

### Key state/data transformations

- Lines 21-23: computes `result` using `result=new_test_result('kernel/perm_group_suite', 'Permutation group metadata', 'real-valued-character Abelian permutation subgroup tables must be exact and internally c…`.
- Lines 26: computes `G` using `G=perm_group('S6A')`.
- Lines 27-28: computes `elements_ref` using `elements_ref=[1 2 3 4 5 6;2 1 4 3 5 6;3 4 1 2 5 6;4 3 2 1 5 6; 1 2 3 4 6 5;2 1 4 3 6 5;3 4 1 2 6 5;4 3 2 1 6 5]`.
- Lines 29-36: computes `chars_ref` using `chars_ref=[1 1 1 1 1 1 1 1; 1 1 -1 -1 1 1 -1 -1; 1 -1 1 -1 1 -1 1 -1; 1 -1 -1 1 1 -1 -1 1; 1 1 1 1 -1 -1 -1 -1; 1 1 -1 -1 -1 -1 1 1; 1 -1 1 -1 -1 1 -1 1; 1 -1 -1 1 -1 1…`.
- Lines 70: computes `group_names` using `group_names={'S4A','S6A','S8A'}`.
- Lines 71: computes `group_degrees` using `group_degrees=[4 6 8]`.
- Lines 72: computes `group_orders` using `group_orders=[4 8 16]`.
- Lines 101: computes `group_closed` using `group_closed=true; group_comm=true; group_invol=true; group_perm=true`.
- Lines 102: computes `identity` using `identity=1:group_degrees(n)`.
- Lines 106: computes `element` using `element=G.elements(k,:)`.
- Lines 107: computes `group_perm` using `group_perm=group_perm&&isequal(sort(element),identity)`.
- Lines 108: computes `group_invol` using `group_invol=group_invol&&isequal(element(element),identity)`.
- Lines 110: computes `prod_left` using `prod_left=G.elements(k,G.elements(m,:))`.
- Lines 111: computes `prod_right` using `prod_right=G.elements(m,G.elements(k,:))`.
- Lines 113: computes `group_comm` using `group_comm=group_comm&&isequal(prod_left,prod_right)`.
- Lines 128: computes `central_count` using `central_count=0; all_perms=perms(identity)`.
- Lines 130: computes `central_perm` using `central_perm=true`.

## Outputs

- result -regression test result with explanatory messages
- The test checks the real-valued-character Abelian permutation subgroups
- against explicit element tables, explicit character tables, character
- orthogonality, closure, commutativity, and involution order.

## Implementation structure

- Tests permutation group database metadata. Syntax:
- result=test_perm_group_suite()
- result -regression test result with explanatory messages
- The test checks the real-valued-character Abelian permutation subgroups
- against explicit element tables, explicit character tables, character
- orthogonality, closure, commutativity, and involution order.
- Announce the test target
- State the utility target of the test
- Check the S6 Abelian subgroup table
- Check the S8 Abelian subgroup table
- Check structural consistency of all real-valued-character Abelian options
- Get the group under test

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `perm_group()`, `test_close()`, `test_true()`, `group_orders()`, `group_degrees()`, `all()`, `isequal()`, `element()`, `ismember()`, `perms()`, `all_perms()`.
