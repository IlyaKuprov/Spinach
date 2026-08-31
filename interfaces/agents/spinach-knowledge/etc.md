# Spinach code index: etc

- Source root: `/home/kuprov/.openclaw/workspace/Spinach`
- Source commit: `f053e432a61d7144f3946d73d0a672e3ccfc3fc5`
- Source tree state: `clean`
- Path set: tracked MATLAB files from `git ls-files '*.m'`; untracked MATLAB files are excluded.
- Files indexed: **47** MATLAB files
- Generated: 2026-08-30T02:52:31

| File | Signature | Summary | LOC |
|---|---|---|---:|
| `etc/data_processing/autophase.m` | `[spec,cheb_coeffs]=autophase(spec,guess)` | Chebyshev phase corrector for 1D NMR spectra. Views the phase profile across the spectral window as a slowly va- rying f | 98 |
| `etc/data_processing/destreak.m` | `spectrum=destreak(spectrum)` | Reduces streak artefacts in 2D and 3D NMR spectra. Edges of the input spectrum must be free of genuine signals. Cell arr | 105 |
| `etc/data_processing/lpredict.m` | `y=lpredict(x,npcoeffs,npredps)` | Forward linear prediction. Syntax: y=lpredict(x,npcoeffs,npredps) Parameters: x -input data, a column vector npcoeffs -n | 79 |
| `etc/diamond_defects/diamond_co.m` | `[sys,inter]=diamond_co(parameters)` | Cobalt-related defect spin system for diamond. Syntax: [sys,inter]=diamond_co(parameters) Magnetic parameters from Nadol | 125 |
| `etc/diamond_defects/diamond_gev0.m` | `[sys,inter]=diamond_gev0(parameters)` | GeV0 spin system for diamond. Syntax: [sys,inter]=diamond_gev0(parameters) Magnetic parameters from Nadolinny et al., Ph | 110 |
| `etc/diamond_defects/diamond_n2vm.m` | `[sys,inter]=diamond_n2vm(parameters)` | N2V-spin system for diamond. Syntax: [sys,inter]=diamond_n2vm(parameters) Magnetic parameters from Green et al., Phys. R | 155 |
| `etc/diamond_defects/diamond_n_inter.m` | `[sys,inter]=diamond_n_inter(parameters)` | Nitrogen interstitial spin system for diamond. Syntax: [sys,inter]=diamond_n_inter(parameters) Magnetic parameters from  | 131 |
| `etc/diamond_defects/diamond_ni.m` | `[sys,inter]=diamond_ni(parameters)` | Nickel-related defect spin system for diamond. Syntax: [sys,inter]=diamond_ni(parameters) W8 magnetic parameters from Lu | 232 |
| `etc/diamond_defects/diamond_nv0_es.m` | `[sys,inter]=diamond_nv0_es(parameters)` | NV0 excited-state spin system for diamond. Syntax: [sys,inter]=diamond_nv0_es(parameters) Magnetic parameters from Felto | 108 |
| `etc/diamond_defects/diamond_nvm_gs.m` | `[sys,inter]=diamond_nvm_gs(parameters)` | NV centre ground state spin system for diamond. Syntax: [sys,inter]=diamond_nvm_gs(parameters) Magnetic parameters from: | 127 |
| `etc/diamond_defects/diamond_ov0.m` | `[sys,inter]=diamond_ov0(parameters)` | Neutral oxygen-vacancy (OV0, WAR5) centre ground state spin system for diamond. Syntax: [sys,inter]=diamond_ov0(paramete | 155 |
| `etc/diamond_defects/diamond_p.m` | `[sys,inter]=diamond_p(parameters)` | Phosphorus-related defect spin system for diamond. Syntax: [sys,inter]=diamond_p(parameters) Magnetic parameters from Na | 170 |
| `etc/diamond_defects/diamond_p1.m` | `[sys,inter]=diamond_p1(parameters)` | P1 centre spin system for diamond. Syntax: [sys,inter]=diamond_p1(parameters) Magnetic parameters from: Nir-Arad et al.  | 116 |
| `etc/diamond_defects/diamond_p1_13c.m` | `[sys,inter]=diamond_p1_13c(parameters)` | P1 centre spin system with 13C neighbours in diamond. Syntax: [sys,inter]=diamond_p1_13c(parameters) The electron and ni | 330 |
| `etc/diamond_defects/diamond_r2.m` | `[sys,inter]=diamond_r2(parameters)` | R2 self-interstitial spin system for diamond. Syntax: [sys,inter]=diamond_r2(parameters) Magnetic parameters from Hunt e | 91 |
| `etc/diamond_defects/diamond_siv0.m` | `[sys,inter]=diamond_siv0(parameters)` | SiV0 spin system for diamond. Syntax: [sys,inter]=diamond_siv0(parameters) Magnetic parameters from Edmonds et al., Phys | 126 |
| `etc/diamond_defects/diamond_ti.m` | `[sys,inter]=diamond_ti(parameters)` | Titanium-related defect spin system for diamond. Syntax: [sys,inter]=diamond_ti(parameters) Magnetic parameters from Nad | 186 |
| `etc/diamond_defects/diamond_vacancy.m` | `[sys,inter]=diamond_vacancy(parameters)` | Vacancy-family defect spin systems for diamond. Syntax: [sys,inter]=diamond_vacancy(parameters) R4/W6 parameters from Tw | 144 |
| `etc/estimators/guess_csa_pro.m` | `CSAs=guess_csa_pro(aa_nums,pdb_ids,coords,options)` | Guesses a reasonable amide bond 15N CSA tensor anisotropy, and a reasonable 13C=O tensor anisotropy, given a local prote | 338 |
| `etc/estimators/guess_j_nuc.m` | `jmatrix=guess_j_nuc(nuc_num,nuc_typ,pdb_id,coords)` | RNA assignments of J-couplings from literature values and Karplus cur- ves. Syntax: jmatrix=guess_j_nuc(nuc_num,nuc_typ, | 462 |
| `etc/estimators/guess_j_pro.m` | `jmatrix=guess_j_pro(aa_num,aa_typ,pdb_id,coords)` | Assigns J-couplings from literature values and Karplus curves. Syntax: jmatrix=guess_j_pro(aa_num,aa_typ,pdb_id,coords)  | 660 |
| `etc/forum.m` | `forum()` | Opens Spinach support forum page. | 18 |
| `etc/hebrew/hebrew.m` | `ncards=hebrew(mode,max_cards) % #NWIKI #NHEAD` | IK's Hebrew flashcards function. The Excel files should contain Hebrew vocabulary in separate spreadsheets: nouns.xlsx - | 389 |
| `etc/mex/compile_mex.m` | `compile_mex() % #NGRUM #NHEAD` | MEX compilation utility. Rebuilds all C++ MEX binaries in the current directory. Syntax: compile_mex() | 42 |
| `etc/molecules/allyl_pyruvate.m` | `[sys,inter]=allyl_pyruvate(spins)` | Spin system of allyl pyruvate. Isotropic chemical shifts and J-couplings determined by spectral fitting, coordinates and | 165 |
| `etc/molecules/cyprinol.m` | `[sys,inter,bas]=cyprinol()` | Spin system of cyprinol. Isotropic chemical shifts and J-couplings are taken from http://dx.doi.org/10.1002/mrc.4782 and | 193 |
| `etc/molecules/dac_reaction.m` | `[sys,inter,bas,kin]=dac_reaction()` | Example Diels-Alder cycloaddition reaction settings: pentadiene (reactant), acrylonitrile (reactant), exo-norbornene (pr | 194 |
| `etc/molecules/fatty_acid.m` | `[sys,inter]=fatty_acid(nprotons)` | Spin system approximating that of a fatty acid. Syntax: [sys,inter]=fatty_acid(nprotons) Parameters: nprotons -the numbe | 71 |
| `etc/molecules/lactate.m` | `[sys,inter]=lactate(spins)` | Spin system of 13C-labelled lactate with the OH protons assumed to be in rapid exchange with water. Syntax: [sys,inter,b | 79 |
| `etc/molecules/methyl_group.m` | `xyz=methyl_group(c_xyz,cc_th,cc_ph,phase)` | Coordinates for the four atoms of a methyl group. Syntax: xyz=methyl_group(c_xyz,cc_th,cc_ph,phase) Parameters: c_xyz -c | 70 |
| `etc/molecules/strychnine.m` | `[sys,inter]=strychnine(spins)` | Spin system of strychnine. Isotropic chemical shifts and J-couplings are taken from "200 and more NMR experiments: a pra | 200 |
| `etc/molecules/zfs_sampling.m` | `[D,E,W]=zfs_sampling(npoints_d,npoints_e,tol)` | Gadolinium ZFS probability distribution function for DOTA-type ligand complexes in cryogenic water-methanol glasses. The | 102 |
| `etc/phantoms/phantoms.m` | `[R1Ph,R2Ph,PDPh,dims,npts]=phantoms(ph_name)` | MRI phantom library. Syntax: [R1Ph,R2Ph,PDPh,dims,npts]=phantoms(ph_name) Parameters: ph_name -character string giving t | 169 |
| `etc/textbook/levelpop.m` | `[E,P,dP]=levelpop(isotope,field,temperature)` | Equilibrium populations of the energy levels of a user-specified spin at the user-specified temperature. Energies are re | 83 |
| `etc/textbook/lorentz.m` | `[J,K,Kil]=lorentz(L)` | The (L,0)(+)(0,L) irreducible matrix representation of the Lorentz group with inversion. Syntax: [J,K,Kil]=lorentz(L) Pa | 83 |
| `etc/textbook/quad_shift.m` | `delta=quad_shift(Cq,eta,v0,S,m)` | Second order shift of the centre of gravity of the powder pattern of \|S,m> to \|S,m-1> transition in the NMR spectrum of  | 73 |
| `etc/textbook/r1csa2tauc.m` | `tauc=r1csa2tauc(R1,del_sq,B0,isotope)` | Estimates the rotational correlation time from the longitudinal CSA relaxation rate. Syntax: tauc=r1csa2tauc(R1,del_sq,B | 67 |
| `etc/textbook/r1n_dnp.m` | `R1n=r1n_dnp(B0,T,g,T1e,T1n_bulk,r,bet)` | A simple model for nuclear longitudnal relaxation rate dies to the presence of an unpaired electron in cryogenic DNP set | 95 |
| `etc/textbook/r2csa2tauc.m` | `tauc=r2csa2tauc(R2,del_sq,B0,isotope)` | Estimates the rotational correlation time from the transverse CSA relaxation rate. Syntax: tauc=r2csa2tauc(R1,del_sq,B0, | 91 |
| `etc/textbook/rlx_csa.m` | `[r1,r2]=rlx_csa(B0,isotope,Z,tau_c)` | Redfield theory expressions for CSA relaxation, including contributions from the antisymmetric part. Syntax: [r1,r2]=rlx | 72 |
| `etc/textbook/rlx_dd_csa.m` | `[A,B,X]=rlx_dd_csa(B0,tau_c,isotopes,deltas,coords)` | Redfield theory expressions for some relaxation and cross- relaxation rates in a CSA-DD-CSA system with two spin-1/2 par | 162 |
| `etc/textbook/rlx_dip.m` | `[r1,r2,rx]=rlx_dip(B0,spins,dist,tau_c)` | Redfield theory expressions for dipolar relaxation and cross- relaxation rates, isotropic tumbling in liquid phase. Synt | 102 |
| `etc/textbook/rlx_hfc.m` | `[r1,r2,rx]=rlx_hfc(B0,HFC,spins,tau_c)` | Redfield theory expressions for hyperfine relaxation and cross- relaxation rates, isotropic tumbling in liquid phase. Sy | 139 |
| `etc/textbook/rlx_nqi.m` | `[r1,r2,t1,t2]=rlx_nqi(I,omega,C_q,eta_q,tau_c)` | Redfield theory expressions for quadrupolar relaxation rates, isotropic tumbling in liquid phase. Syntax: [r1,r2,t1,t2]= | 88 |
| `etc/textbook/rlx_sbm.m` | `[r1,r2]=rlx_sbm(B0,nucleus,dist,a_iso,e_spin,g_eff,t1e,t2e,tau_r)` | Solomon-Bloembergen-Morgan nuclear relaxation rates due to a paramagnetic centre. Syntax: [r1,r2]=rlx_sbm(B0,nucleus,dis | 131 |
| `etc/textbook/trosy_eff.m` | `eff=trosy_eff(B0,isotopes,xyz,csa)` | TROSY efficiency in a two-spin system. Returns the extent of the cancellation of the CSA contribution to the trans- vers | 111 |
| `etc/wiki.m` | `wiki()` | Opens Spinach documentation Wiki page. | 21 |
