# interfaces/gaussian/gslice.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/interfaces/gaussian/gslice.m`
- Signature: `gslice()`
- Total lines: 103

## Purpose

Slices a Gaussian geometry scan log into property calculation inputs at the energy minimum geometries. The function asks for the log and for a text file containing the property calcula- tion header. Some paths are hard-coded; edit as appropriate.

## Physical / mathematical content

- Gaussian interfaces. These parse quantum-chemistry output into spin Hamiltonian ingredients such as hyperfine, shielding, or exchange parameters.

## Numerical / algorithmic content

- The file also defines local helper function(s): `dump()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 12-19: Assign atomic symbols; implemented by `periodic_table={'H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg','Al','Si','P','S','Cl','Ar', 'K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','A…`.
- Lines 21-22: Read the log file; implemented by `[g03_output,g03_output_path]=uigetfile('*.*','Relaxed geometry scan log')`.
- Lines 25-26: Locate and the stationary point reports; implemented by `minima=[]`.
- Lines 34-35: Locate standard orientation entry points; implemented by `std_geom_start=[]`.
- Lines 46-47: Locate standard orientation end points; implemented by `std_geom_end=[]`.
- Lines 58-59: Read the standard orientations; implemented by `coord_blocks=cell(length(std_geom_start),1)`.
- Lines 67-68: Get the header; implemented by `[g03_header,g03_header_path]=uigetfile('*.*','Gaussian header')`.
- Lines 71-72: Write the inputs; implemented by `runscript=fopen('compute.bat','a')`.

### Control flow inferred from the code

- Line 27: `for` loop over `n=1:length(g03_output)`.
- Line 28: conditional branch on `strcmp(deblank(g03_output(n)),'-- Stationary point found.')`.
- Line 36: `for` loop over `n=minima`.
- Line 37: `for` loop over `k=n:-1:1`.
- Line 38: conditional branch on `strcmp(deblank(g03_output(k)),'Standard orientation:')`.
- Line 48: `for` loop over `n=std_geom_start`.
- Line 49: `for` loop over `k=n:length(g03_output)`.
- Line 50: conditional branch on `strcmp(deblank(g03_output(k)),'---------------------------------------------------------------------')`.
- Line 60: `for` loop over `n=1:length(std_geom_start)`.
- Line 62: `for` loop over `k=std_geom_start(n):std_geom_end(n)`.
- Line 76: `for` loop over `n=1:length(coord_blocks)`.
- Line 79: `for` loop over `k=1:size(coord_blocks{n},1)`.

### Key state/data transformations

- Lines 13-19: computes `periodic_table` using `periodic_table={'H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg','Al','Si','P','S','Cl','Ar', 'K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','A…`.
- Lines 22: computes `[g03_output,g03_output_path]` using `[g03_output,g03_output_path]=uigetfile('*.*','Relaxed geometry scan log')`.
- Lines 23: computes `g03_output` using `g03_output=textread([g03_output_path g03_output],'%s','bufsize',65536,'delimiter','\n')`.
- Lines 26: computes `minima` using `minima=[]`.
- Lines 35: computes `std_geom_start` using `std_geom_start=[]`.
- Lines 47: computes `std_geom_end` using `std_geom_end=[]`.
- Lines 59: computes `coord_blocks` using `coord_blocks=cell(length(std_geom_start),1)`.
- Lines 61: computes `coord_blocks{n}` using `coord_blocks{n}=zeros(std_geom_end(n)-std_geom_start(n),6)`.
- Lines 63: computes `coord_blocks{n}(k-std_geom_start(n)+1,:)` using `coord_blocks{n}(k-std_geom_start(n)+1,:)=cell2mat(textscan(g03_output{k},'%n %n %n %n %n %n'))`.
- Lines 68: computes `[g03_header,g03_header_path]` using `[g03_header,g03_header_path]=uigetfile('*.*','Gaussian header')`.
- Lines 69: computes `header` using `header=textread([g03_header_path g03_header],'%s','delimiter','\n')`.
- Lines 72: computes `runscript` using `runscript=fopen('compute.bat','a')`.
- Lines 74: computes `fprintf(runscript,'%s\n','set GAUSS_EXEDIR` using `fprintf(runscript,'%s\n','set GAUSS_EXEDIR=C:\G16W')`.
- Lines 75: computes `fprintf(runscript,'%s\n','set GAUSS_SCRDIR` using `fprintf(runscript,'%s\n','set GAUSS_SCRDIR=C:\Temp')`.
- Lines 77: computes `g03_input` using `g03_input=fopen(['g16_input_' num2str(n) '.gjf'],'a')`.
- Lines 80: computes `current_line` using `current_line=[periodic_table{coord_blocks{n}(k,2)} ' ' num2str(coord_blocks{n}(k,4:6),'%1.8f ')]`.

### Local helper functions

- Line 92: `dump()` — `function dump(file_id,cell_array)`. Don't show a half-done job to a fool. A Russian saying
  - Representative operation: `for p=1:length(cell_array)`.
  - Representative operation: `fprintf(file_id,'%s\n',cell_array{p})`.

## Implementation structure

- Slices a Gaussian geometry scan log into property calculation
- inputs at the energy minimum geometries. The function asks for
- the log and for a text file containing the property calcula-
- tion header. Some paths are hard-coded; edit as appropriate.
- Assign atomic symbols
- Read the log file
- Locate and the stationary point reports
- Locate standard orientation entry points
- Locate standard orientation end points
- Read the standard orientations
- Get the header
- Write the inputs

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `uigetfile()`, `textread()`, `strcmp()`, `deblank()`, `g03_output()`, `num2str()`, `std_geom_end()`, `std_geom_start()`, `cell2mat()`, `textscan()`, `fopen()`, `dump()`, `fclose()`.
