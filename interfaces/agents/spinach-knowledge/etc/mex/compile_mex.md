# etc/mex/compile_mex.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/etc/mex/compile_mex.m`
- Signature: `compile_mex() % #NGRUM #NHEAD`
- Total lines: 42

## Purpose

MEX compilation utility. Rebuilds all C++ MEX binaries in the current directory. Syntax: compile_mex()

## Physical / mathematical content

- This file belongs to the `etc` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content


## Code-derived implementation details

### Comment-guided execution stages

- Lines 12-13: Get own location; implemented by `P=mfilename('fullpath'); P=P(1:(end-20))`.
- Lines 15-17: Lorentzian convolution; implemented by `mex('-R2018a','-O','-DNDBUG','COMPFLAGS=$COMPFLAGS','LINKFLAGS=$LINKFLAGS', [P '/kernel/line_shapes/lorentzcon.cpp'],'-outdir',[P '/kernel/line_shapes'])`.
- Lines 19-21: Gaussian convolution; implemented by `mex('-R2018a','-O','-DNDBUG','COMPFLAGS=$COMPFLAGS','LINKFLAGS=$LINKFLAGS', [P '/kernel/line_shapes/gausscon.cpp'],'-outdir',[P '/kernel/line_shapes'])`.
- Lines 23-25: Cubic polynomial roots; implemented by `mex('-R2018a','-O','-DNDBUG','COMPFLAGS=$COMPFLAGS','LINKFLAGS=$LINKFLAGS', [P '/kernel/eigenfields/cubic_roots.cpp'],'-outdir',[P '/kernel/eigenfields'])`.
- Lines 27-29: Sparse double row sorter; implemented by `mex('-R2018a','-O','-DNDEBUG', [P '/kernel/indexing/spsortrows.cpp'],'-outdir',[P '/kernel/indexing'])`.
- Lines 31-33: Sparse double unique columns; implemented by `mex('-R2018a','-O','-DNDEBUG', [P '/kernel/indexing/spunicols.cpp'],'-outdir',[P '/kernel/indexing'])`.

### Key state/data transformations

- Lines 13: computes `P` using `P=mfilename('fullpath'); P=P(1:(end-20))`.
- Lines 16-17: computes `mex('-R2018a','-O','-DNDBUG','COMPFLAGS` using `mex('-R2018a','-O','-DNDBUG','COMPFLAGS=$COMPFLAGS','LINKFLAGS=$LINKFLAGS', [P '/kernel/line_shapes/lorentzcon.cpp'],'-outdir',[P '/kernel/line_shapes'])`.

## Implementation structure

- MEX compilation utility. Rebuilds all C++ MEX binaries in
- the current directory. Syntax:
- compile_mex()
- Get own location
- Lorentzian convolution
- Gaussian convolution
- Cubic polynomial roots
- Sparse double row sorter
- Sparse double unique columns
- Audiophiles don't use their equip-
- ment to listen to your music. Audio-
- philes use your music to listen to

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `mfilename()`, `mex()`.
