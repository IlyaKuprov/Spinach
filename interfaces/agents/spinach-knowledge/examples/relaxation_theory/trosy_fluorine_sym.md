# examples/relaxation_theory/trosy_fluorine_sym.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/relaxation_theory/trosy_fluorine_sym.m`
- Signature: `trosy_fluorine_sym()`
- Total lines: 74

## Purpose

Transverse relaxation rate as a function of the applied magnetic field in a 3-fluorotyrosine labelled protein. The fluorine atom and its directly bonded carbon are included. Analytical calcula- tions broken down by mechanism. Calculation time: seconds.

## Physical / mathematical content

- Relaxation-theory examples. The mathematical backbone is Bloch-Redfield-Wangsness or stochastic Liouville theory, spectral densities, cross-correlation terms, motional models, and extraction of longitudinal/transverse decay behaviour from superoperators.
- Chemical-shift anisotropy is present: shielding is treated as a second-rank tensor whose orientation relative to the field or rotor axis modulates line shapes and transfer dynamics.

## Numerical / algorithmic content


## Code-derived implementation details

### Comment-guided execution stages

- Lines 12-14: Read 3-fluorotyrosine DFT calculation; implemented by `[~,inter_dft]=g2spinach(gparse('../standard_systems/3_fluoro_tyr.log'), {{'C','13C'},{'F','19F'}},[186.38 192.97])`.
- Lines 16-17: Extract coordinates and CSAs; implemented by `ZF=inter_dft.zeeman.matrix{8}`.
- Lines 22-23: Magnetic field grid; implemented by `lin_freq=linspace(200,800,20)`.
- Lines 26-27: Loop over magnetic fields; implemented by `for n=1:numel(B0)`.
- Lines 29-30: Call the analytical function; implemented by `[F,C]=rlx_dd_csa(B0(n),25e-9,{'19F','13C'},{ZF,ZC},{RF,RC})`.
- Lines 32-33: Relaxation rates; implemented by `r2c(n)=C.r2.total`.
- Lines 40-41: Mechanisms for 13C; implemented by `c_tro_dd(n)=C.trosy.dd`.
- Lines 47-48: Plotting; implemented by `kfigure()`.
- Lines 64-65: TROSY rate by mechanism; implemented by `kfigure()`.

### Control flow inferred from the code

- Line 27: `for` loop over `n=1:numel(B0)`.

### Key state/data transformations

- Lines 13-14: computes `[~,inter_dft]` using `[~,inter_dft]=g2spinach(gparse('../standard_systems/3_fluoro_tyr.log'), {{'C','13C'},{'F','19F'}},[186.38 192.97])`.
- Lines 17: computes `ZF` using `ZF=inter_dft.zeeman.matrix{8}`.
- Lines 18: computes `ZC` using `ZC=inter_dft.zeeman.matrix{7}`.
- Lines 19: computes `RF` using `RF=inter_dft.coordinates{8}`.
- Lines 20: computes `RC` using `RC=inter_dft.coordinates{7}`.
- Lines 23: computes `lin_freq` using `lin_freq=linspace(200,800,20)`.
- Lines 24: computes `B0` using `B0=2*pi*lin_freq*1e6/spin('1H')`.
- Lines 30: computes `[F,C]` using `[F,C]=rlx_dd_csa(B0(n),25e-9,{'19F','13C'},{ZF,ZC},{RF,RC})`.
- Lines 33: computes `r2c(n)` using `r2c(n)=C.r2.total`.
- Lines 34: computes `r2f(n)` using `r2f(n)=F.r2.total`.
- Lines 35: computes `f_bro(n)` using `f_bro(n)=F.trosy.total_bro`.
- Lines 36: computes `f_nar(n)` using `f_nar(n)=F.trosy.total_nar`.
- Lines 37: computes `c_bro(n)` using `c_bro(n)=C.trosy.total_bro`.
- Lines 38: computes `c_nar(n)` using `c_nar(n)=C.trosy.total_nar`.
- Lines 41: computes `c_tro_dd(n)` using `c_tro_dd(n)=C.trosy.dd`.
- Lines 42: computes `c_tro_csa(n)` using `c_tro_csa(n)=C.trosy.csa`.
- Lines 43: computes `c_tro_xc(n)` using `c_tro_xc(n)=C.trosy.xc`.

## Implementation structure

- Transverse relaxation rate as a function of the applied magnetic
- field in a 3-fluorotyrosine labelled protein. The fluorine atom
- and its directly bonded carbon are included. Analytical calcula-
- tions broken down by mechanism.
- Calculation time: seconds.
- Read 3-fluorotyrosine DFT calculation
- Extract coordinates and CSAs
- Magnetic field grid
- Loop over magnetic fields
- Call the analytical function
- Relaxation rates
- Mechanisms for 13C

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `g2spinach()`, `gparse()`, `spin()`, `rlx_dd_csa()`, `r2c()`, `r2f()`, `f_bro()`, `f_nar()`, `c_bro()`, `c_nar()`, `c_tro_dd()`, `c_tro_csa()`, `c_tro_xc()`, `kfigure()`, `ylim()`, `kxlabel()`.
