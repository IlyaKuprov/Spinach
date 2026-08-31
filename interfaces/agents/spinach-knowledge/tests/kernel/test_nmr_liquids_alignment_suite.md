# tests/kernel/test_nmr_liquids_alignment_suite.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_nmr_liquids_alignment_suite.m`
- Signature: `result=test_nmr_liquids_alignment_suite()`
- Total lines: 193

## Purpose

Tests compact literature-alignment probes for liquid-state NMR pulse sequences. Syntax: result=test_nmr_liquids_alignment_suite()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 17-18: Announce the test target; implemented by `fprintf('TESTING: NMR liquids literature-alignment probes\n')`.
- Lines 20-23: State the test target; implemented by `result=new_test_result('kernel/nmr_liquids_alignment_suite', 'NMR liquids literature-alignment probes', 'Updated liquid-state pulse sequence paths must produce finite no…`.
- Lines 25-26: Build a compact homonuclear system for gCOSY; implemented by `sys.magnet=9.4`.
- Lines 35-36: Set compact gCOSY parameters; implemented by `parameters.sweep=100`.
- Lines 46-47: Run the P-type gCOSY pathway; implemented by `fid_p=liquid(spin_system,@gcosy,parameters,'nmr')`.
- Lines 49-50: Run the N-type gCOSY pathway; implemented by `parameters.pathway='N'`.
- Lines 53-54: Run the phase-sensitive gCOSY pathway pair; implemented by `parameters.pathway='P+N'`.
- Lines 57-59: Check gCOSY outputs; implemented by `result=test_close(result,'gCOSY P finite',all(isfinite(fid_p(:))),true,0,0, 'the P-type gradient pathway must produce finite data')`.
- Lines 77-78: Recombine the gCOSY pathways for phase-sensitive processing; implemented by `f1_pos=fftshift(fft(fid_pn.pos,parameters.npoints(2),1),1)`.
- Lines 83-86: Check the recombined gCOSY spectrum; implemented by `result=test_close(result,'gCOSY P+N recombined finite', all(isfinite(spec_pn(:)))&&(norm(real(spec_pn(:)))>0),true,0,0, 'P+N pathway recombination must produce finite no…`.
- Lines 88-89: Build a compact heteronuclear system for HMBC and HSQC; implemented by `sys.isotopes={'1H','15N'}`.
- Lines 95-96: Set compact HMBC parameters; implemented by `parameters=struct()`.
- Lines 103-104: Run HMBC on a non-carbon heteronucleus; implemented by `fid_hmbc=liquid(spin_system,@hmbc,parameters,'nmr')`.
- Lines 106-109: Check HMBC output; implemented by `result=test_close(result,'HMBC non-carbon finite', all(isfinite(fid_hmbc(:)))&&(norm(fid_hmbc(:))>0),true,0,0, 'HMBC coherence selection must follow parameters.spins{1}')`.
- Lines 111-112: Set compact HSQC parameters; implemented by `parameters.decouple_f1={'1H'}`.
- Lines 115-116: Run HSQC; implemented by `fid_hsqc=liquid(spin_system,@hsqc,parameters,'nmr')`.
- Lines 118-121: Check HSQC outputs; implemented by `result=test_close(result,'HSQC positive pathway finite', all(isfinite(fid_hsqc.pos(:)))&&(norm(fid_hsqc.pos(:))>0),true,0,0, 'HSQC positive pathway must produce finite o…`.
- Lines 126-127: Run CT-HSQC; implemented by `fid_ct_hsqc=liquid(spin_system,@ct_hsqc,parameters,'nmr')`.

### Key state/data transformations

- Lines 21-23: computes `result` using `result=new_test_result('kernel/nmr_liquids_alignment_suite', 'NMR liquids literature-alignment probes', 'Updated liquid-state pulse sequence paths must produce finite no…`.
- Lines 26: computes `sys.magnet` using `sys.magnet=9.4`.
- Lines 27: computes `sys.isotopes` using `sys.isotopes={'1H','1H'}`.
- Lines 28: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={0 0}`.
- Lines 29: computes `inter.coupling.scalar` using `inter.coupling.scalar=cell(2,2)`.
- Lines 30: computes `inter.coupling.scalar{1,2}` using `inter.coupling.scalar{1,2}=8`.
- Lines 31: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 32: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 33: computes `spin_system` using `spin_system=test_spin_system(sys,inter,bas)`.
- Lines 36: computes `parameters.sweep` using `parameters.sweep=100`.
- Lines 37: computes `parameters.npoints` using `parameters.npoints=[4 4]`.
- Lines 38: computes `parameters.spins` using `parameters.spins={'1H'}`.
- Lines 39: computes `parameters.angle` using `parameters.angle=pi/2`.
- Lines 40: computes `parameters.g_amp` using `parameters.g_amp=3`.
- Lines 41: computes `parameters.g_dur` using `parameters.g_dur=1e-3`.
- Lines 42: computes `parameters.g_stab_del` using `parameters.g_stab_del=0`.
- Lines 43: computes `parameters.s_len` using `parameters.s_len=1.5`.
- Lines 44: computes `parameters.pathway` using `parameters.pathway='P'`.

## Outputs

- result -regression test result with explanatory messages
- The test covers compact gCOSY, HMBC, HSQC, and TOCSY paths that
- were updated during the nmr_liquids literature-alignment pass.

## Implementation structure

- Tests compact literature-alignment probes for liquid-state NMR
- pulse sequences. Syntax:
- result=test_nmr_liquids_alignment_suite()
- result -regression test result with explanatory messages
- The test covers compact gCOSY, HMBC, HSQC, and TOCSY paths that
- were updated during the nmr_liquids literature-alignment pass.
- Announce the test target
- State the test target
- Build a compact homonuclear system for gCOSY
- Set compact gCOSY parameters
- Run the P-type gCOSY pathway
- Run the N-type gCOSY pathway

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `test_spin_system()`, `liquid()`, `test_close()`, `all()`, `fid_p()`, `fid_n()`, `isstruct()`, `isfield()`, `fftshift()`, `conj()`, `spec_pn()`, `fid_hmbc()`, `isequal()`, `state()`, `speye()`.
