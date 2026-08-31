# tests/interfaces/test_orca_parser.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/interfaces/test_orca_parser.m`
- Signature: `result=test_orca_parser()`
- Total lines: 53

## Purpose

Tests the ORCA log parser on the logs bundled with the examples. Syntax: result=test_orca_parser()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content


## Code-derived implementation details

### Comment-guided execution stages

- Lines 19-20: Announce the test target; implemented by `fprintf('TESTING: ORCA log parser\n')`.
- Lines 22-25: State the physical target of the test; implemented by `result=new_test_result('interfaces/orca_parser', 'ORCA log parser', 'magnetic parameters are read from ORCA logs and attached to the correct atoms.')`.
- Lines 27-28: Locate the example logs bundled with Spinach; implemented by `root_dir=fileparts(fileparts(fileparts(mfilename('fullpath'))))`.
- Lines 30-32: Methyl radical, a vacuum DFT calculation with a g-tensor and hyperfines; implemented by `props=oparse(fullfile(root_dir,'examples','esr_liq_pulsed', 'data_import','orca_methyl_radical.out'))`.
- Lines 43-44: Copper porphyrin, hyperfine tensors printed for the protons only; implemented by `props=oparse(fullfile(root_dir,'examples','visualisation','porphyrine.out'))`.

### Key state/data transformations

- Lines 23-25: computes `result` using `result=new_test_result('interfaces/orca_parser', 'ORCA log parser', 'magnetic parameters are read from ORCA logs and attached to the correct atoms.')`.
- Lines 28: computes `root_dir` using `root_dir=fileparts(fileparts(fileparts(mfilename('fullpath'))))`.
- Lines 31-32: computes `props` using `props=oparse(fullfile(root_dir,'examples','esr_liq_pulsed', 'data_import','orca_methyl_radical.out'))`.
- Lines 38: computes `proton_hfc` using `proton_hfc=gauss2mhz(props.hfc.full.matrix{2})`.
- Lines 45: computes `tensor_atoms` using `tensor_atoms=find(~cellfun(@isempty,props.hfc.full.matrix))`.
- Lines 46: computes `proton_atoms` using `proton_atoms=find(strcmp(props.symbols,'H'))`.

## Outputs

- result -regression test result with explanatory messages
- The test checks version detection, the mapping of ORCA's zero-based
- and incompletely printed nucleus numbering onto the atoms of the
- coordinate table, and the numerical values of the g-tensor and of
- the hyperfine tensors against the quantities that ORCA itself
- prints separately in the same logs.

## Implementation structure

- Tests the ORCA log parser on the logs bundled with the examples. Syntax:
- result=test_orca_parser()
- result -regression test result with explanatory messages
- The test checks version detection, the mapping of ORCA's zero-based
- and incompletely printed nucleus numbering onto the atoms of the
- coordinate table, and the numerical values of the g-tensor and of
- the hyperfine tensors against the quantities that ORCA itself
- prints separately in the same logs.
- Announce the test target
- State the physical target of the test
- Locate the example logs bundled with Spinach
- Methyl radical, a vacuum DFT calculation with a g-tensor and hyperfines

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `fileparts()`, `mfilename()`, `oparse()`, `fullfile()`, `test_true()`, `strcmp()`, `test_close()`, `gauss2mhz()`, `cellfun()`, `isequal()`.
