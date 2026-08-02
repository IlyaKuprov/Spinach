# Literature

Theory implemented in Spinach, and representative published work that used it.
Format is `Authors, "Title", Journal Volume, Pages/Article (Year), DOI`; an
absent field was not available from a checkable record and must not be filled
in from memory.

## Citing the library

The library paper is the citation for Spinach itself; cite it together with the
methodology papers for the machinery that carried the calculation - basis
restriction, Fokker-Planck spatial dynamics, relaxation theory, optimal control.

- H.J. Hogben, M. Krzystyniak, G.T.P. Charnock, P.J. Hore, I. Kuprov, "Spinach - a software library for simulation of spin dynamics in large spin systems", Journal of Magnetic Resonance 208(2), 179-194 (2011), DOI: 10.1016/j.jmr.2010.11.008
- I. Kuprov, "Large-scale NMR simulations in liquid state: a tutorial", Magnetic Resonance in Chemistry 56(6), 415-437 (2018), DOI: 10.1002/mrc.4660
- I. Kuprov, "Defeating the Matrix", Journal of Magnetic Resonance 306, 75-79 (2019), DOI: 10.1016/j.jmr.2019.07.031

## State space restriction and basis construction

The reduced Liouville-space basis is what makes large spin systems tractable.

- I. Kuprov, N. Wagner-Rundell, P.J. Hore, "Polynomially scaling spin dynamics simulation algorithm based on adaptive state-space restriction", Journal of Magnetic Resonance 189(2), 241-250 (2007), DOI: 10.1016/j.jmr.2007.09.014 - founding paper behind `sphten-liouv` and the `IK-*` bases.
- I. Kuprov, "Polynomially scaling spin dynamics II: further state-space compression using Krylov subspace techniques and zero track elimination", Journal of Magnetic Resonance 195(1), 45-51 (2008), DOI: 10.1016/j.jmr.2008.08.008 - zero track elimination, in `zte.m`.
- H.J. Hogben, P.J. Hore, I. Kuprov, "Strategies for state space restriction in densely coupled spin systems with applications to spin chemistry", Journal of Chemical Physics 132(17), 174101 (2010), DOI: 10.1063/1.3398146 - symmetry factorisation in Liouville space, in `symmetry.m`.
- A. Karabanov, I. Kuprov, G.T.P. Charnock, A. van der Drift, L.J. Edwards, W. Kockenberger, "On the accuracy of the state space restriction approximation for spin dynamics simulations", Journal of Chemical Physics 135(8), 084106 (2011), DOI: 10.1063/1.3624564 - error bounds and basis selection criteria, cited in `basis.m`.
- M. Krzystyniak, L.J. Edwards, I. Kuprov, "Destination state screening of active spaces in spin dynamics simulations", Journal of Magnetic Resonance 210(2), 228-232 (2011), DOI: 10.1016/j.jmr.2011.03.010 - backward reachability screening, in `path_trace.m`.

## Large matrix representations

Tensor-structured and unopened-Kronecker forms for matrices too large to store.

- D.V. Savostyanov, S.V. Dolgov, J.M. Werner, I. Kuprov, "Exact NMR simulation of protein-size spin systems using tensor train formalism", Physical Review B 90(8), 085139 (2014), DOI: 10.1103/PhysRevB.90.085139 - the `@ttclass` overload.
- L.J. Edwards, D.V. Savostyanov, Z.T. Welderufael, D. Lee, I. Kuprov, "Quantum mechanical NMR simulation algorithm for protein-size spin systems", Journal of Magnetic Resonance 243, 107-113 (2014), DOI: 10.1016/j.jmr.2014.04.002 - restricted basis plus tensor trains at protein scale.
- A.J. Allami, M.G. Concilio, P. Lally, I. Kuprov, "Quantum mechanical MRI simulations: solving the matrix dimension problem", Science Advances 5(7), eaaw8962 (2019), DOI: 10.1126/sciadv.aaw8962 - the `@polyadic` and `@opium` classes.
- I.V. Oseledets, "Tensor-Train Decomposition", SIAM Journal on Scientific Computing 33(5), 2295-2317 (2011), DOI: 10.1137/090752286 - the underlying decomposition, cited in `ttclass.m`.

## Fokker-Planck formalism and spatial dynamics

Spinning, diffusion, flow and gradients enter as a spin-space direct product.

- I. Kuprov, "Fokker-Planck formalism in magnetic resonance simulations", Journal of Magnetic Resonance 270, 124-135 (2016), DOI: 10.1016/j.jmr.2016.07.005 - the generator structure behind the MAS, diffusion, flow and imaging contexts.
- A. Karabanov, A. van der Drift, L.J. Edwards, I. Kuprov, W. Kockenberger, "Quantum mechanical simulation of solid effect dynamic nuclear polarisation using Krylov-Bogolyubov time averaging and a restricted state-space", Physical Chemistry Chemical Physics 14(8), 2658-2668 (2012), DOI: 10.1039/C2CP23233B - doubly-rotating-frame time averaging, in `average.m`.

## Propagation

- L.J. Edwards, I. Kuprov, "Parallel density matrix propagation in spin dynamics simulations", Journal of Chemical Physics 136(4), 044108 (2012), DOI: 10.1063/1.3679656 - parallel propagation of state-vector stacks, in `evolution.m`.
- D.L. Goodwin, I. Kuprov, "Auxiliary matrix formalism for interaction representation transformations, optimal control, and spin relaxation theories", Journal of Chemical Physics 143(8), 084113 (2015), DOI: 10.1063/1.4928978 - propagator caching, rotating-frame transformations, directional derivatives, Redfield integrals.

## Relaxation theory

Bloch-Redfield-Wangsness superoperators built without diagonalising the
Hamiltonian, and the correlation functions that feed them.

- I. Kuprov, "Diagonalization-free implementation of spin relaxation theory for large spin systems", Journal of Magnetic Resonance 209(1), 31-38 (2011), DOI: 10.1016/j.jmr.2010.12.004 - the algorithm in `redfield_integral_{serial,async}.m`.
- I. Kuprov, L.C. Morris, J.N. Glushka, J.H. Prestegard, "Using molecular dynamics trajectories to predict nuclear spin relaxation behaviour in large spin systems", Journal of Magnetic Resonance 323, 106891 (2021), DOI: 10.1016/j.jmr.2020.106891 - relaxation superoperators from MD trajectories.
- M.G. Concilio, M. Soundararajan, L. Frydman, I. Kuprov, "High-field solution state DNP using cross-correlations", Journal of Magnetic Resonance 326, 106940 (2021), DOI: 10.1016/j.jmr.2021.106940 - the dipolar/CSA cross-correlation path.
- M.G. Concilio, I. Kuprov, L. Frydman, "J-driven dynamic nuclear polarization for sensitizing high field solution state NMR", Physical Chemistry Chemical Physics 24(4), 2118-2125 (2022), DOI: 10.1039/d1cp04186j - scalar-coupling-driven cross-relaxation DNP.

## Optimal control

GRAPE-family pulse optimisation: gradients, Hessians, waveforms, distortions.

- I. Kuprov, C.T. Rodgers, "Derivatives of spin dynamics simulations", Journal of Chemical Physics 131(23), 234108 (2009), DOI: 10.1063/1.3267086 - propagator-derivative theory, superseded by auxiliary matrices.
- P. de Fouquieres, S.G. Schirmer, S.J. Glaser, I. Kuprov, "Second order gradient ascent pulse engineering", Journal of Magnetic Resonance 212(2), 412-417 (2011), DOI: 10.1016/j.jmr.2011.07.023 - the BFGS/LBFGS GRAPE optimisers.
- I. Kuprov, "Spin system trajectory analysis under optimal control pulses", Journal of Magnetic Resonance 233, 107-112 (2013), DOI: 10.1016/j.jmr.2013.02.012 - `trajan.m` and `trajsimil.m`.
- D.L. Goodwin, I. Kuprov, "Modified Newton-Raphson GRAPE methods for optimal control of spin systems", Journal of Chemical Physics 144(20), 204107 (2016), DOI: 10.1063/1.4949534 - the Hessian-regularised optimiser, `fmaxnewton.m`.
- M.S. Vinding, D.L. Goodwin, I. Kuprov, T.E. Lund, "Optimal control gradient precision trade-offs: application to fast generation of DeepControl libraries for MRI", Journal of Magnetic Resonance 333, 107094 (2021), DOI: 10.1016/j.jmr.2021.107094 - reduced-precision gradients for large pulse libraries.
- U. Rasulov, A. Acharya, M. Carravetta, G. Mathies, I. Kuprov, "Simulation and design of shaped pulses beyond the piecewise-constant approximation", Journal of Magnetic Resonance 353, 107478 (2023), DOI: 10.1016/j.jmr.2023.107478 - trapezium waveform machinery, `trapdiff.m`.
- U. Rasulov, I. Kuprov, "Instrumental distortions in quantum optimal control", Journal of Chemical Physics 162(16), 164107 (2025), DOI: 10.1063/5.0264092 - amplifier and filter models in `optimcon/distortions/`.
- I. Kuprov, "Optimal control of spin systems", in Spin: From Basic Symmetries to Quantum Optimal Control, Springer, 313-349 (2023), DOI: 10.1007/978-3-031-05607-9_8 - book-chapter treatment of the stack.
- O.W. Sorensen, "Polarization transfer experiments in high-resolution NMR spectroscopy", Progress in Nuclear Magnetic Resonance Spectroscopy 21(6), 503-569 (1989), DOI: 10.1016/0079-6565(89)80006-8 - universal bound on transferable polarisation, in `sorensen.m`.

## Fitting and inverse problems

- E.A. Suturina, D. Haussinger, K. Zimmermann, L. Garbuio, M. Yulikov, G. Jeschke, I. Kuprov, "Model-free extraction of spin label position distributions from pseudocontact shift data", Chemical Science 8(4), 2751-2757 (2017), DOI: 10.1039/c6sc03736d - Tikhonov-regularised reconstruction with L-curve parameter selection.

## Applications: hyperpolarisation, PHIP and SABRE

- J. Eronen et al., "Characterisation of the polarisation transfer to fluorinated pyridines in SABRE", Physical Chemistry Chemical Physics (2025), DOI: 10.1039/d5cp01418b
- A. Ortmeier et al., "SABRE-hyperpolarization dynamics of [1-13C]pyruvate monitored by in situ zero- to ultra-low field NMR", Journal of Magnetic Resonance Open (2024), DOI: 10.1016/j.jmro.2024.100149
- S. Lehmkuhl et al., "SABRE polarized low field rare-spin spectroscopy", Journal of Chemical Physics (2020), DOI: 10.1063/5.0002412

## Applications: dynamic nuclear polarisation

- S.A. Jegadeesan et al., "Simulation of pulsed dynamic nuclear polarization in the steady state", Journal of Chemical Physics (2025), DOI: 10.1063/5.0283196
- V.S. Redrouthu et al., "Efficient Pulsed Dynamic Nuclear Polarization with the X-Inverse-X Sequence", Journal of the American Chemical Society (2022), DOI: 10.1021/jacs.1c09900
- F.A. Perras et al., "Full-Scale Ab Initio Simulation of Magic-Angle-Spinning Dynamic Nuclear Polarization", Journal of Physical Chemistry Letters (2020), DOI: 10.1021/acs.jpclett.0c00955

## Applications: solid-state NMR and magic angle spinning

- S. Ray et al., "Optimal control-based nuclear spin cross-polarization in the presence of complicating anisotropic interactions", Physical Chemistry Chemical Physics (2025), DOI: 10.1039/d5cp00096c
- E. Burlinson et al., "Significant acceleration of solid-state NMR simulations via three-angle powder averaging", Journal of Magnetic Resonance (2025), DOI: 10.1016/j.jmr.2025.107966
- B. Simoes de Almeida et al., "Theory and simulations of homonuclear three-spin systems in rotating solids", Journal of Chemical Physics (2021), DOI: 10.1063/5.0055583

## Applications: EPR, DEER and spin labels

- J. Keeley et al., "Neural networks in pulsed dipolar spectroscopy: a practical guide", Journal of Magnetic Resonance (2022), DOI: 10.1016/j.jmr.2022.107186
- J.D. Lehner et al., "Modeling of motional EPR spectra using hindered Brownian rotational diffusion and the stochastic Liouville equation", Journal of Chemical Physics (2020), DOI: 10.1063/1.5139935
- N. Manukovsky et al., "Time domain simulation of Gd3+-Gd3+ distance measurements by EPR", Journal of Chemical Physics (2017), DOI: 10.1063/1.4994084

## Applications: paramagnetic NMR

- T. Muntener et al., "Pseudocontact Shifts in Biomolecular NMR Spectroscopy", Chemical Reviews (2022), DOI: 10.1021/acs.chemrev.1c00796
- H.W. Orton et al., "Accurate Electron-Nucleus Distances from Paramagnetic Relaxation Enhancements", Journal of the American Chemical Society (2018), DOI: 10.1021/jacs.8b03858
- J. Rantaharju et al., "Spin dynamics simulation of electron spin relaxation in Ni2+(aq)", Journal of Chemical Physics (2014), DOI: 10.1063/1.4885050

## Applications: radical pairs and spin chemistry

- I.V. Zhukov et al., "Simulation of electron and nuclear spin dynamics in many-spin charge-separated states", Journal of Chemical Physics (2025), DOI: 10.1063/5.0244106
- L.M. Antill et al., "RadicalPy: A Tool for Spin Dynamics Simulations", Journal of Chemical Theory and Computation (2024), DOI: 10.1021/acs.jctc.4c00887
- C. Nielsen et al., "MolSpin - flexible and extensible general spin dynamics software", Journal of Chemical Physics (2019), DOI: 10.1063/1.5125043

## Applications: singlet order and long-lived states

- B.B. Kharkov et al., "Weak nuclear spin singlet relaxation mechanisms revealed by experiment and computation", Physical Chemistry Chemical Physics (2022), DOI: 10.1039/d1cp05537b
- Y. Feng et al., "Long-lived polarization protected by symmetry", Journal of Chemical Physics (2014), DOI: 10.1063/1.4896895
- Y. Feng et al., "Accessing long-lived nuclear singlet states between chemically equivalent spins without breaking symmetry", Nature Physics (2012), DOI: 10.1038/nphys2425

## Applications: zero, ultralow and Earth field NMR

- J.E. Elenewski et al., "Prospects for NMR spectral prediction on fault-tolerant quantum computers", Physical Review Research (2026), DOI: 10.1103/km66-ltpl
- I. Mandzhieva et al., "Zero-Field NMR and Millitesla-SLIC Spectra for >200 Molecules from Density Functional Theory and Spin Dynamics", Journal of Chemical Information and Modeling (2025), DOI: 10.1021/acs.jcim.5c00111
- D.C. Kaseman et al., "Earth's Field NMR for Organophosphate Chemical Warfare Agent Detection", Applied Magnetic Resonance (2023), DOI: 10.1007/s00723-023-01565-4

## Applications: imaging, flow and diffusion

- R. Mishra et al., "Theoretical analysis of flow effects in spatially encoded diffusion NMR", Journal of Chemical Physics (2022), DOI: 10.1063/5.0130125
- K. Landheer et al., "Theoretical description of modern 1H in vivo magnetic resonance spectroscopic pulse sequences", Journal of Magnetic Resonance Imaging (2019), DOI: 10.1002/jmri.26846
- E.W. Zhao et al., "In situ NMR metrology reveals reaction mechanisms in redox flow batteries", Nature (2020), DOI: 10.1038/s41586-020-2081-7

## Applications: pulse design and quantum simulation

- M. Foroozandeh, "Spin dynamics during chirped pulses: applications to homonuclear decoupling and broadband excitation", Journal of Magnetic Resonance (2020), DOI: 10.1016/j.jmr.2020.106768
- J. Saywell et al., "Biselective pulses for large-area atom interferometry", Physical Review A (2020), DOI: 10.1103/physreva.101.063625
- M.G. Algaba et al., "Co-Design quantum simulation of nanoscale NMR", Physical Review Research (2022), DOI: 10.1103/physrevresearch.4.043089
