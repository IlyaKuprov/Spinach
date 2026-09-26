# Spinach code index: kernel

- Source repository: `https://github.com/IlyaKuprov/Spinach`

| File | Signature | Summary |
|---|---|---|
| `kernel/assume.m` | `spin_system=assume(spin_system,assumptions,retention)` | Sets case-specific assumptions for various simulation contexts. This function determines the behaviour of the Hamiltonia |
| `kernel/average.m` | `H=average(spin_system,Hp,H0,Hm,omega,theory)` | Average Hamiltonian theories under Zeeman interaction rotating frame transformations. Syntax: H=average(spin_system,Hp,H |
| `kernel/basis.m` | `spin_system=basis(spin_system,bas)` | Basis set control. This is the second mandatory function (after create.m) that must be called in every calculation to bu |
| `kernel/cache/bos_product_table.m` | `[product_table_left,product_table_right]=bos_product_table(nlevels)` | Structure coefficient tables for the associative envelopes of truncated Weyl algebras spanned by orthogonalised bosonic |
| `kernel/cache/cacheman.m` | `cacheman(spin_system) %#NHEAD` | Cache management heuristics. Looks after the scratch folder and prevents it from filling up the disk. Do not call direct |
| `kernel/cache/ist_product_table.m` | `[product_table_left,product_table_right]=ist_product_table(mult)` | Structure coefficient tables for the associative envelopes of su(mult) algebras. Syntax: [product_table_left,product_tab |
| `kernel/cache/sle_operators.m` | `[Lx,Ly,Lz,D,space_basis]=sle_operators(max_rank,int_ranks)` | Wigner D function basis set and rotation generators required by the SLE module. Syntax: [Lx,Ly,Lz,D,space_basis]=sle_ope |
| `kernel/cache/st_product_table.m` | `[pt_left,pt_right]=st_product_table(nlevels)` | Structure coefficient tables for single transition operators. Syntax: [pt_left,pt_right]=st_product_table(nlevels) Param |
| `kernel/cache/wipe_cache.m` | `wipe_cache(spin_system)` | Forces a wipe of the Spinach cache folder. Syntax: wipe_cache(spin_system) Parameters: spin_system -Spinach object with |
| `kernel/carrier.m` | `H=carrier(spin_system,spins,operator_type)` | Returns the "carrier" Hamiltonian -the part of the Zeeman interaction Hamiltonian that corresponds to all particles havi |
| `kernel/coherence.m` | `rho=coherence(spin_system,rho,spec)` | Coherence order selection function -keeps only the specified orders of coherence in the state vector. This is useful as |
| `kernel/coherent.m` | `rho=coherent(spin_system,mode,alpha)` | Coherent state of a bosonic mode. Builds the normalised trunca- tion of the coherent state with the specified amplitude |
| `kernel/contexts/crystal.m` | `answer=crystal(spin_system,pulse_sequence,parameters,assumptions)` | Single-crystal interface to pulse sequences. Generates a Liouvillian superoperator and passes it on to the pulse sequenc |
| `kernel/contexts/device.m` | `answer=device(spin_system,pulse_sequence,parameters,assumptions)` | Spin-boson device interface to pulse sequences. Generates the evolution generators for a device containing spins and bos |
| `kernel/contexts/doublerot.m` | `[answer,sph_grid]=doublerot(spin_system,pulse_sequence,parameters,assumptions)` | Double angle spinning context. In Liouville space, this wrapper builds the Fokker-Planck evolution generator and passes |
| `kernel/contexts/floquet.m` | `[answer,sph_grid]=floquet(spin_system,pulse_sequence,parameters,assumptions)` | Floquet magic angle spinning context. Generates a Liouvillian super- operator and passes it on to the pulse sequence fun |
| `kernel/contexts/gridfree.m` | `answer=gridfree(spin_system,pulse_sequence,parameters,assumptions)` | Fokker-Planck magic angle spinning and SLE context. Generates a Liouvil- lian superoperator and passes it on to the puls |
| `kernel/contexts/imaging.m` | `answer=imaging(spin_system,pulse_sequence,parameters)` | Fokker-Planck imaging simulation context. Generates the Hamiltonian, the relaxation superoperator, the kinetics superope |
| `kernel/contexts/liquid.m` | `answer=liquid(spin_system,pulse_sequence,parameters,assumptions)` | Liquid-phase interface to pulse sequences. Generates a Liouvillian superoperator and passes it on to the pulse sequence |
| `kernel/contexts/meshflow.m` | `answer=meshflow(spin_system,pulse_sequence,parameters)` | First draft of the magnetohydrodynamics context for microfluidic simu- lations. Generates evolution generators and passe |
| `kernel/contexts/powder.m` | `[answer,sph_grid]=powder(spin_system,pulse_sequence,parameters,assumptions)` | Static powder interface to pulse sequences. Generates a Liouvillian superoperator, the initial state, the coil state, th |
| `kernel/contexts/singlerot.m` | `[answer,sph_grid]=singlerot(spin_system,pulse_sequence,parameters,assumptions)` | Single angle spinning context. In Liouville space, this wrapper builds the Fokker-Planck evolution generator that includ |
| `kernel/conventions/transforms/anas2mat.m` | `M=anas2mat(iso,an,as,alp,bet,gam)` | Converts anisotropy and asymmetry representation of a 3x3 interaction tensor (Haeberlen-Mehring convention) into the cor |
| `kernel/conventions/transforms/anax2dcm.m` | `dcm=anax2dcm(rot_axis,rot_angle)` | Converts angle-axis rotation parameters to a direction cosine matrix in the active convention, matching the one used by |
| `kernel/conventions/transforms/anax2qter.m` | `q=anax2qter(rot_axis,rot_angle)` | Converts angle-axis rotation parameters into a quaternion. Syntax: q=anax2qter(rot_axis,rot_angle) Parameters: rot_axis |
| `kernel/conventions/transforms/ang2cgsppm.m` | `cgsppm=ang2cgsppm(ang)` | Converts magnetic susceptibility from the Angstrom^3 units required by Spinach pseudocontact shift functionality into th |
| `kernel/conventions/transforms/axis_tsymm.m` | `A=axis_tsymm(T,a)` | Roughly averages an interaction tensor with respect to the rotation around a user-specified axis. Syntax: T=axis_tsymm(T |
| `kernel/conventions/transforms/axrh2mat.m` | `M=axrh2mat(iso,ax,rh,alp,bet,gam)` | Converts axiality and rhombicity representation of a 3x3 interaction tensor into the corresponding matrix. Syntax: M=axr |
| `kernel/conventions/transforms/cart2mode.m` | `mode_derivs=cart2mode(cart_derivs,eigvecs,masses,frqs)` | Converts Cartesian derivatives of spin Hamiltonian parameters, as produced by electronic structure theory packages, into |
| `kernel/conventions/transforms/castep2nqi.m` | `nqi=castep2nqi(V,Q,I)` | Converts CASTEP EFG tensor (it is printed in atomic units) to NQI 3x3 tensor in Hz that is required by Spinach. Syntax: |
| `kernel/conventions/transforms/cgsppm2ang.m` | `ang=cgsppm2ang(cgsppm)` | Converts magnetic susceptibility from the cgs-ppm (aka cm^3/mol) units quoted by quantum chemistry packages into Angstro |
| `kernel/conventions/transforms/dcm2euler.m` | `[arg1,arg2,arg3]=dcm2euler(dcm)` | Converts directional cosine matrix into Euler angles, ZYZ active convention (rotating the object rather than the axes). |
| `kernel/conventions/transforms/dcm2qter.m` | `q=dcm2qter(dcm)` | Converts a direction cosine matrix in the active convention of euler2dcm.m function into a unit quaternion. Syntax: q=dc |
| `kernel/conventions/transforms/dcm2wigner.m` | `D=dcm2wigner(dcm)` | Converts a directional cosine matrix into second-rank Wigner function matrix. Syntax: D=dcm2wigner(dcm) Parameters: dcm |
| `kernel/conventions/transforms/eeqq2nqi.m` | `Q=eeqq2nqi(C_q,eta_q,I,eulers)` | Converts the C_q and eta_q quadrupolar interaction specification convention into a 3x3 interaction matrix in Hz. Syntax: |
| `kernel/conventions/transforms/ejec2duffing.m` | `[frq,anharm]=ejec2duffing(ej,ec)` | Converts the Josephson and charging energies of a transmon into the Duffing oscillator frequency and anharmonicity expec |
| `kernel/conventions/transforms/euler2dcm.m` | `R=euler2dcm(arg1,arg2,arg3)` | Converts Euler angles (ZYZ active convention) into a direction cosine matrix. Syntax: R=euler2dcm(alpha,beta,gamma) OR R |
| `kernel/conventions/transforms/euler2qter.m` | `q=euler2qter(arg1,arg2,arg3)` | Converts Euler angles (ZYZ active convention) into a unit quaternion in the active convention, matching euler2dcm.m func |
| `kernel/conventions/transforms/euler_equiv.m` | `answer=euler_equiv(eulers_a,eulers_b,tol)` | Checks whether two ZYZ active Euler angle sets specify the same rotation. Syntax: answer=euler_equiv(eulers_a,eulers_b,t |
| `kernel/conventions/transforms/euler_sup.m` | `rot_cmp=euler_sup(rot_one,rot_two)` | Superposition of ZYZ active Euler rotations. Syntax: rot_cmp=euler_sup(rot_one,rot_two) Parameters: rot_one -first Euler |
| `kernel/conventions/transforms/frac2cart.m` | `[XYZ,va,vb,vc]=frac2cart(a,b,c,alp,bet,gam,ABC)` | Converts fractional crystallographic coordinates to Cartesian coordinates. Syntax: [XYZ,va,vb,vc]=frac2cart(a,b,c,alpha, |
| `kernel/conventions/transforms/fwhm2rlx.m` | `r2rate=fwhm2rlx(fwhm)` | Converts full width at half-maximum (FWHM) of an NMR signal into an approximation of the R2 rate. Syntax: r2rate=fwhm2rl |
| `kernel/conventions/transforms/g2freq.m` | `f=g2freq(g,B)` | Converts g-tensor units into electron Zeeman frequency units. Syntax: f=g2freq(g,B) Parameters: g - g-values, scalar or |
| `kernel/conventions/transforms/gauss2mhz.m` | `hfc_mhz=gauss2mhz(hfc_gauss,g)` | Converts hyperfine couplings from Gauss to MHz (linear frequency). The Gauss specification may be defined as "the magnet |
| `kernel/conventions/transforms/ham2nqi.m` | `[omega,Q]=ham2nqi(H)` | Converts a single-spin Hamiltonian back into the Zeeman and quadrupolar interaction parameters that had been used to gen |
| `kernel/conventions/transforms/hartree2joule.m` | `energy=hartree2joule(energy)` | Converts Hartree energy units into J/mol. A Hartree is twice the ground state ionisation energy of the hydrogen atom. Sy |
| `kernel/conventions/transforms/hz2icm.m` | `icm=hz2icm(hz)` | Converts Hz units used in magnetic resonance into cm^-1 units used in spectroscopy. Syntax: icm=hz2icm(hz) Arrays of any |
| `kernel/conventions/transforms/hz2ppm.m` | `ppm=hz2ppm(hz,B0,nucleus)` | Converts resonance offsets into chemical shifts. Syntax: ppm=hz2ppm(hz,B0,nucleus) Parameters: hz -resonance offset in H |
| `kernel/conventions/transforms/ias2mat.m` | `C=ias2mat(a,d,A)` | Reconstruction of a 3x3 real interaction matrix C between real vectors u and v from its isotropic-antisymmetric-symmetri |
| `kernel/conventions/transforms/icm2hz.m` | `hz=icm2hz(icm)` | Converts cm^-1 units used in spectroscopy into Hz units preferred in magnetic resonance. Syntax: hz=icm2hz(icm) Arrays o |
| `kernel/conventions/transforms/kelvin2hz.m` | `hz=kelvin2hz(kelvin)` | Converts Kelvin energy units used for Debye temperatures and thermal energy scales in solid state physics into Hz units |
| `kernel/conventions/transforms/mat2axrh.m` | `[iso,ax,rh,eigvals]=mat2axrh(M)` | Computes axiality and rhombicity of a symmetric 3x3 interaction tensor from the corresponding matrix. Syntax: [iso,ax,rh |
| `kernel/conventions/transforms/mat2ias.m` | `[a,d,A]=mat2ias(C)` | Isotropic-antisymmetric-symmetric decomposition of a 3x3 real interaction matrix between real vectors u and v: u'*C*v = |
| `kernel/conventions/transforms/mat2sphten.m` | `[rank0,rank1,rank2]=mat2sphten(M)` | Converts a 3x3 interaction matrix into the irreducible spherical tensor notation: one rank 0 component, three rank 1 com |
| `kernel/conventions/transforms/mev2hz.m` | `hz=mev2hz(mev)` | Converts meV energy units used in solid state physics and phonon spectroscopy into Hz units preferred in magnetic resona |
| `kernel/conventions/transforms/mhz2gauss.m` | `hfc_gauss=mhz2gauss(hfc_mhz,g)` | Converts hyperfine couplings from MHz (linear frequency) to Gauss. The Gauss specification may be defined as "the magnet |
| `kernel/conventions/transforms/mt2hz.m` | `hfc_hz=mt2hz(hfc_mt,g)` | Converts hyperfine couplings from milliTesla to Hz (linear frequency). The milliTesla specification may be defined as "t |
| `kernel/conventions/transforms/ppm2hz.m` | `hz=ppm2hz(ppm,B0,nucleus)` | Converts chemical shifts into resonance offsets. Syntax: hz=ppm2hz(ppm,B0,nucleus) Parameters: ppm -chemical shift in pp |
| `kernel/conventions/transforms/qform2sph.m` | `[r0,r1,r2]=qform2sph(A)` | Returns the spherical harmonic expansion coefficients of the following quadratic form: [x y z]*A*[x y z]'/norm([x y z],2 |
| `kernel/conventions/transforms/qter2anax.m` | `[rot_axis,rot_angle]=qter2anax(q)` | Converts a quaternion representation of a rotation into angle-axis rotation parameters. Syntax: [rot_axis,rot_angle]=qte |
| `kernel/conventions/transforms/qter2dcm.m` | `dcm=qter2dcm(q)` | Converts a unit quaternion into a direction cosine matrix in the active convention, matching the one used by euler2dcm.m |
| `kernel/conventions/transforms/qter2euler.m` | `[alpha,beta,gamma]=qter2euler(q)` | Converts a unit quaternion in the active convention into Euler angles (ZYZ active convention), matching euler2dcm.m func |
| `kernel/conventions/transforms/rotmat_align.m` | `rot_mat=rotmat_align(v_from,v_to)` | Rotation matrix aligning one vector with another vector. Syntax: rot_mat=rotmat_align(v_from,v_to) Parameters: v_from -t |
| `kernel/conventions/transforms/sphten2mat.m` | `M=sphten2mat(rank0,rank1,rank2)` | Converts the nine components of the irreducible spherical tensor re- presentation of an interaction tensor into the Cart |
| `kernel/conventions/transforms/spsk2mat.m` | `M=spsk2mat(iso,sp,sk,alp,bet,gam)` | Converts span and skew representation of a 3x3 interaction tensor (Herzfeld-Berger convention) into the corresponding ma |
| `kernel/conventions/transforms/stev2sph.m` | `Bkq=stev2sph(k,Bkq)` | Transforms the coefficients in front of Stevens operators, as produced by stevens.m, into the coefficients before the ir |
| `kernel/conventions/transforms/tsm2param.m` | `[ax,rh,angles]=tsm2param(M)` | Attempts to convert a traceless symmetric 3x3 interaction matrix into axiality, rhombicity and three Euler angles. The t |
| `kernel/conventions/transforms/weblab2nqi.m` | `varargout=weblab2nqi(C_q,eta_q,I,alpha,theta,phi)` | Converts the Weblab one-cone model parameters (see weblab_cone.png) into NQI tensors used by Spinach. Syntax: [Q1,Q2]=we |
| `kernel/conventions/transforms/xyz2sph.m` | `[r,theta,phi] = xyz2sph(x,y,z)` | Converts Cartesian coordinates [x y z] into spherical coordinates according to the ISO convention. Syntax: [r,theta,phi] |
| `kernel/conventions/transforms/zfs2mat.m` | `M=zfs2mat(D,E,alp,bet,gam)` | Converts D and E zero-field splitting parameters described in the abstract of (http://dx.doi.org/10.1063/1.1682294) into |
| `kernel/correlation.m` | `rho=correlation(spin_system,rho,orders,spins)` | Correlation order selection function -keeps only the specified orders of spin correlation in the state vector. This is u |
| `kernel/create.m` | `spin_system=create(sys,inter)` | The entry function of the Spinach kernel that creates the spin system object that the rest of the library requires to ru |
| `kernel/decouple.m` | `[L,rho]=decouple(spin_system,L,rho,spins)` | Obliterates all interactions and populations in the subspace of states that involve the specified spins in any way. The |
| `kernel/derivatives/fdhess.m` | `H=fdhess(A,nstenc)` | Returns the finite-difference Hessian of a 3D array using a finite difference scheme with a user-specified number of ste |
| `kernel/derivatives/fdkup.m` | `K=fdkup(npoints,extents,chi,nstenc)` | Returns a finite difference representation of the Kuprov operator: K[rho]=-(1/3)*Trace(Hessian[rho]*chi) with the number |
| `kernel/derivatives/fdlap.m` | `L=fdlap(dims,extents,nstenc)` | Returns a finite-difference representation of the Laplacian for an array with a user-specified finite difference stencil |
| `kernel/derivatives/fdmat.m` | `D=fdmat(dim,nstenc,order,boundary)` | Returns arbitrary-order central finite-difference differentiation matrices (sparse) with unit grid point spacing. Syntax |
| `kernel/derivatives/fdvec.m` | `dx=fdvec(x,npoints,order)` | Performs arbitrary-order finite-difference differentiation of a user-supplied row or column vector. Uses central finite- |
| `kernel/derivatives/fdweights.m` | `w=fdweights(target_point,grid_points,max_order)` | Calculates finite difference weights for numerical derivatives, including order 0, which amounts to interpolation. Synta |
| `kernel/derivatives/fftdiff.m` | `kern=fftdiff(order,npoints,dx)` | Spectral differentiation kernel. Syntax: kern=fftdiff(order,npoints,dx) Parameters: order -order of the derivative npoin |
| `kernel/derivatives/fourdif.m` | `[x,DM]=fourdif(N,m)` | The function [x,DM] = fourdif(N,m) computes the m'th derivative Fourier spectral differentiation matrix on grid with N e |
| `kernel/derivatives/fourlap.m` | `L=fourlap(npoints,extents)` | Returns a Fourier spectral representation of the Laplacian acting on a 3D data array. Syntax: L=fourlap(npoints,extents) |
| `kernel/derivatives/pseudomodulation.m` | `output=pseudomodulation(field,spectrum,mod_amp,mod_order)` | Pseudomodulation of uniformly sampled spectra using the Hyde et al. Fourier-domain algorithm. Syntax: output=pseudomodul |
| `kernel/derivatives/sgolaydiff.m` | `dy=sgolaydiff(y,der_order,npoints,poly_order)` | Savitzky-Golay differentiation of noisy sampled signals by local least-squares polynomial fitting. Syntax: dy=sgolaydiff |
| `kernel/eigenfields/cubic_roots.m` | `root_list=cubic_roots(poly_coeffs,root_tol)` | Real roots of a cubic polynomial in the unit interval. Syntax: root_list=cubic_roots(poly_coeffs,root_tol) Parameters: p |
| `kernel/eigenfields/eigenfields.m` | `tran=eigenfields(spin_system,parameters,Hz,Hc,Hmw)` | Computes resonance fields. For a Hamiltonian Hc+B*Hz, returns all magnetic fields B for which the difference between two |
| `kernel/eigenfields/rootmatch.m` | `[idx1,idx2,idx3]=rootmatch(field1,field2,field3,edge12,edge23,edge31)` | Global order-preserving root matching between three magnetic field root lists. The routine returns indices into the inpu |
| `kernel/eigenfields/voitlander.m` | `spec=voitlander(spin_system,parameters,triangle,Ic,Iz,Qc,Qz,Hmw)` | Adaptively recursed Voitlander integrator. Computes an approximation of an integral of field-swept EPR transition over a |
| `kernel/evolution.m` | `answer=evolution(spin_system,L,coil,rho,timestep,nsteps,output,destination)` | Time evolution function. Performs all types of time propagation with automatic trajectory level state space restriction. |
| `kernel/frqoffset.m` | `H=frqoffset(spin_system,H,parameters)` | Adds omega*Lz Larmor frequency offsets to the Hamiltonian; this is useful in liquid state NMR experiments. Syntax: H=frq |
| `kernel/grids/arclength.m` | `sig=arclength(r1,r2)` | Arc length between two points on the unit sphere specified by the unit vectors supplied. Syntax: sig=arclength(r1,r2) Pa |
| `kernel/grids/gaussleg.m` | `[x,w]=gaussleg(a,b,n)` | Computes Gauss-Legendre points and weights in [a,b] interval with accuracy order n. Syntax: [x,w]=gaussleg(a,b,n) Parame |
| `kernel/grids/get_hull.m` | `[hull,edges]=get_hull(theta_angles,phi_angles)` | Generates a convex hull of a two-angle grid for 2D surface plotting. Syntax: [hull,edges]=get_hull(theta_angles,phi_angl |
| `kernel/grids/grid_fibon.m` | `[alps,bets,gams,whts,vorn]=grid_fibon(type,parm)` | Fibonacci type spherical quadrature grids, as per Appendix A.5 of http://dx.doi.org/10.1016/j.jmr.2014.05.009 Syntax: [a |
| `kernel/grids/grid_igloo.m` | `[alps,bets,gams,whts,vorn]=grid_igloo(n_long)` | Igloo grid, as per Appendix A2 of Syntax: [alps,bets,gams,whts,vorn]=grid_igloo(n_long) Parameters: n_long -number of lo |
| `kernel/grids/grid_kron.m` | `[angles,weights]=grid_kron(angles1,weights1,angles2,weights2)` | Spherical grid direct product. Tiles one grid using the rotations of the other. Grids should be supplied using Euler ang |
| `kernel/grids/grid_plot.m` | `grid_plot(x,y,z,vorn,c,options)` | Spherical quadrature grid plotter. Takes a cloud of points on a sphere and plots its Voronoi tessellation. Syntax: grid_ |
| `kernel/grids/grid_polar.m` | `[phi,r,L]=grid_polar(ncircles,rmax)` | Generates a balanced polar grid in which the density of points does not increase towards the centre. Syntax: [phi,r,L]=g |
| `kernel/grids/grid_test.m` | `grid_profile=grid_test(alphas,betas,gammas,weights,ranks,sfun)` | Plots grid integration quality as a function of spherical rank. The quality is defined as the norm of the residual of sp |
| `kernel/grids/grid_trian.m` | `[alps,bets,gams,whts,vorn]=grid_trian(type,n)` | Triangular spherical quadrature grids, as per Appendix A.6 of (http://dx.doi.org/10.1016/j.jmr.2014.05.009). Syntax: [al |
| `kernel/grids/ngridpts.m` | `n=ngridpts(grad_amps,grad_durs,isotope,max_coh_order,sample_size)` | Estimates the minimum number of spatial grid points necessary to have a valid treatment of gradient driven experiments w |
| `kernel/grids/one_vcell_solidangle.m` | `S=one_vcell_solidangle(v,centre)` | Solid angle of a convex spherical polygon as described in Syntax: A=one_vcell_solidangle(v,centre) Parameters: v -(3 x n |
| `kernel/grids/repulsion.m` | `[alphas,betas,gammas,weights]=repulsion(npoints,ndims,niter)` | Generates repulsion grids on a unit hypersphere. See the paper by Bak and Nielsen (http://dx.doi.org/10.1006/jmre.1996.1 |
| `kernel/grids/shrewd.m` | `weights=shrewd(alphas,betas,gammas,max_rank,max_error)` | Computes SHREWD weights for a given two-or three-angle spherical grid. See the paper by Eden and Levitt for details on n |
| `kernel/grids/sphtarea.m` | `S=sphtarea(r1,r2,r3,sflag)` | Area of the curvilinear triangle on the unit sphere defined by the vertex coordinates supplied. Syntax: S=sphtarea(r1,r2 |
| `kernel/grids/sphtrsubd.m` | `[r12,r23,r31]=sphtrsubd(r1,r2,r3)` | Spherical triangle subdivision. Returns the midpoints of the sides of a spherical triangle specified by the unit vectors |
| `kernel/grids/vcell_solidangle.m` | `S=vcell_solidangle(P,K,xyz)` | Solid angle of a spherical Voronoi cell. Syntax: s=vcell_solidangle(P,K,xyz) Parameters: P -(3 x m) array with coordinat |
| `kernel/grids/voronoisphere.m` | `[vertices,indices,polygons,sangles]=voronoisphere(xyz)` | Voronoi tessellation of the unit sphere around the specified po- ints, computed as the exact geometric dual of the Delau |
| `kernel/hamiltonian.m` | `[I,Q]=hamiltonian(spin_system,operator_type)` | Hamiltonian operator or superoperator and its rotational decomposi- tion. Descriptor and operator generation are paralle |
| `kernel/homospoil.m` | `rho=homospoil(spin_system,rho,zqc_flag)` | Emulates a strong homospoil pulse -only zero-frequency states with respect to the carrier frequencies (chemical shifts a |
| `kernel/includes/autoexec.m` | `(script file)` | This include is executed at the start of create.m, it over- rides all user input. A good use case is forcing polyadic or |
| `kernel/includes/end_disallow_gpu.m` | `(script file)` | Reinstates GPU arithmetic setting to the previous state after the start_disallow_gpu command had been issued. |
| `kernel/includes/parallel_profiler_report.m` | `(script file)` | An include that writes the report of the profiling infrastructure around parallel stages. Should be invoked just after a |
| `kernel/includes/parallel_profiler_start.m` | `(script file)` | An include that starts profiling infrastructure around parallel stages. Should be invoked just before a parfor or an spm |
| `kernel/includes/redfield_integral_async.m` | `job_id=brw_compute_kernel(spin_system,w,job_id,upper_lim)` | Bloch-Wangsness-Redfield and Nakajima-Zwanzig integral evaluation, the asynchronous parallel path. This include is calle |
| `kernel/includes/redfield_integral_serial.m` | `(script file)` | Bloch-Wangsness-Redfield and Nakajima-Zwanzig integral evaluati- on, the serial path. This include is called from within |
| `kernel/includes/start_disallow_gpu.m` | `(script file)` | Forces GPU arithmetic to be turned off even if the user had requested it in sys.enable setting; restore previous state u |
| `kernel/indexing/kq2lin.m` | `I=kq2lin(N,K,Q,idx_base)` | Converts k,q indexing of matrices into their linear serpentine indexing. In base 1 indexing convention: (1,1)(1,2)(1,3) |
| `kernel/indexing/lin2kq.m` | `[K,Q]=lin2kq(N,I,idx_base)` | Converts linear serpentine indexing of matrices into their k,q indexing. In base-1 indexing convention: (1,1)(1,2)(1,3) |
| `kernel/indexing/lin2lm.m` | `[L,M]=lin2lm(I)` | Converts linear indexing of spin states into L,M indexing. In the linear indexing convention, spin states are listed in |
| `kernel/indexing/lin2lmn.m` | `[L,M,N]=lin2lmn(I)` | Converts linear indices of Wigner D functions into L,M,N indices. In the linear indexing convention, Wigner D functions |
| `kernel/indexing/lm2lin.m` | `I=lm2lin(L,M)` | Converts L,M indexing of spin states into linear indexing. In the linear indexing convention, spin states are listed in |
| `kernel/indexing/lmn2lin.m` | `I=lmn2lin(L,M,N)` | Converts L,M,N indices of Wigner D functions into linear indices. In the linear indexing convention, Wigner D functions |
| `kernel/indexing/serpentine.m` | `S=serpentine(nlevels,idx_base)` | Serpentine index matrix used in Spinach for single-index numbering of matrix elements. Syntax: S=serpentine(nlevels,idx_ |
| `kernel/indexing/spsortrows.m` | `idx=spsortrows(A)` | Sparse matrix row-sorting permutation utility. Syntax: idx=spsortrows(A) Parameters: A -sparse real double matrix Output |
| `kernel/indexing/spunicols.m` | `A=spunicols(A)` | Sparse matrix unique-column utility. Syntax: A=spunicols(A) Parameters: A -sparse real double matrix Outputs: A -sparse |
| `kernel/integrity/existentials.m` | `existentials()` | Kernel integrity control. Checks for collisions between Spinach functions and anything else that the user may have insta |
| `kernel/integrity/exorcise.m` | `exorcise(mode)` | Searches Spinach distribution folders for any functions that do not conform to the house style. Opens the first one and |
| `kernel/integrity/patrol.m` | `patrol(test_subject)` | This function runs contiouously on one of our servers, its purpose is to catch any unintended consequences before they p |
| `kernel/integrity/rearm.m` | `rearm()` | Rearms the sniffer database. The sniffer checks Spinach distribution .m files for any modifications that the user did si |
| `kernel/integrity/smack.m` | `smack()` | Gives Matlab a good smack every time MDCS gets its kni- ckers in a twist. Syntax: smack() This function shuts down the p |
| `kernel/integrity/sniff.m` | `sniff(action)` | Kernel integrity control. Checks Spinach distribution .m files for any modifications that the user did since downloading |
| `kernel/kinetics/equilibrate.m` | `c=equilibrate(K,c0)` | Equilibrates linear chemical kinetics and returns a vector of equilibrium concentrations. Syntax: c=equilibrate(K,c0) Pa |
| `kernel/kinetics/flow_gen.m` | `F=flow_gen(spin_system,parameters)` | Hydrodynamic flow generator on a mesh. Builds diffusion and flow generator using the mesh parameters in the spin_system |
| `kernel/kinetics/kinetics.m` | `K=kinetics(spin_system)` | Chemical kinetics superoperator. Syntax: K=kinetics(spin_system) Parameters: spin_system - Spinach spin system descripti |
| `kernel/kinetics/react_gen.m` | `G=react_gen(spin_system,reaction)` | Chemical reaction generator builder. Syntax: G=react_gen(spin_system,reaction) Parameters: reaction.reactants -a vector |
| `kernel/krylov.m` | `answer=krylov(spin_system,L,coil,rho,timestep,nsteps,output)` | Krylov propagation function. Avoids matrix exponentiation, but can be slow. Should be used when the Liouvillian exponent |
| `kernel/line_shapes/dhofun.m` | `y=dhofun(x,nat_freq,fwhm)` | Normalised damped harmonic oscillator response function in mag- netic resonance notation. This is the standard shape of |
| `kernel/line_shapes/gausscon.m` | `y=gausscon(offs,ampl,fwhm,x)` | Normalised Gaussian function in magnetic resonance notation and its convolution with a triangular function. Syntax: y=ga |
| `kernel/line_shapes/gaussfun.m` | `y=gaussfun(x,fwhm)` | Normalized Gaussian function in magnetic resonance notation. Syntax: y=gaussfun(x,fwhm) Parameters: x -argument values, |
| `kernel/line_shapes/lorentzcon.m` | `y=lorentzcon(offs,ampl,fwhm,x)` | Normalised Lorentzian function in magnetic resonance notation and its convolution with a triangular function. Syntax: y= |
| `kernel/line_shapes/lorentzfun.m` | `[real_part,imag_part]=lorentzfun(offs,ampl,fwhm,x,phi)` | Normalised Lorentzian function in magnetic resonance notation with a phase distortion. Syntax: [real_part,imag_part]=lor |
| `kernel/multiprop.m` | `rho=multiprop(spin_system,P,rho,N)` | Applies a propagator repeatedly by binary adaptive squaring. Syntax: rho=multiprop(spin_system,P,rho,N) Parameters: spin |
| `kernel/operator.m` | `A=operator(spin_system,operators,spins,operator_type,format)` | Generates Hilbert space operators or Liouville space superoperators from their human-readable descriptions. Syntax: A=op |
| `kernel/operators/bos2ist.m` | `[states,coeffs]=bos2ist(prod_spec,nlevels)` | Irreducible spherical tensor expansion of a user-specified bosonic operator product. Syntax: [states,coeffs]=bos2ist(pro |
| `kernel/operators/boson_mono.m` | `B=boson_mono(nlevels)` | Bosonic monomial operators of the following structure: B(k,q)=(Cr^k)*(An^q) obeying the following commutation relations |
| `kernel/operators/boson_ortho.m` | `B=boson_ortho(nlevels)` | Orthogonal bosonic monomials calculated from the bosonic mono- mial basis produced by boson_mono(nlevels). Gram-Schmidt |
| `kernel/operators/centrans.m` | `A=centrans(mult,type)` | Central transition operators of half-integer spins in the Pauli basis. Syntax: A=centrans(mult,type) Parameters: mult -m |
| `kernel/operators/ct2ist.m` | `[states,coeffs]=ct2ist(mult,type)` | Irreducible spherical tensor expansion of central transition operators of half-integer spins. Syntax: [states,coeffs]=ct |
| `kernel/operators/enlev2bm.m` | `[states,coeffs]=enlev2bm(nlevels,lvl_num)` | Bosonic monomial expansion of specific bosonic energy level projectors. Syntax: [states,coeffs]=enlev2bm(nlevels,lvl_num |
| `kernel/operators/enlev2ist.m` | `[states,coeffs]=enlev2ist(mult,lvl_num,particle)` | Irreducible spherical tensor expansion of specific Zeeman energy level projectors. Syntax: [states,coeffs]=enlev2ist(mul |
| `kernel/operators/hilb2liouv.m` | `L=hilb2liouv(H,conv_type)` | Converts Hilbert space operators into Liouville space super- operators or state vectors. Syntax: L=hilb2liouv(H,conv_typ |
| `kernel/operators/irr_sph_ten.m` | `T=irr_sph_ten(mult,k)` | Single-spin irreducible spherical tensor operators T(k,m) obeying the following commutation relation: [Lz,T_km]=m*T_km S |
| `kernel/operators/lindbladian.m` | `R=lindbladian(A_left,A_right,rho,rlx_rate)` | Generates a Lindblad superoperator from user-specified left-side and right-side product superoperators and calibrates it |
| `kernel/operators/mprealloc.m` | `A=mprealloc(spin_system,nnzpc)` | Preallocates an operator in the current basis. Syntax: A=mprealloc(spin_system,nnzpc) Parameters: nnzpc - expected numbe |
| `kernel/operators/oper2bm.m` | `[states,coeffs]=oper2bm(A)` | Bosonic monomial operator expansion of a user-specified square matrix. Syntax: [states,coeffs]=oper2bm(A) Parameters: A |
| `kernel/operators/oper2ist.m` | `[states,coeffs]=oper2ist(A)` | Irreducible spherical tensor operator expansion of a user- specified square matrix. Syntax: [states,coeffs]=oper2ist(A) |
| `kernel/operators/pauli.m` | `S=pauli(mult)` | Pauli spin operators (sparse, see below for normalisa- tion conventions) for a spin of a user-specified ener- gy level m |
| `kernel/operators/sin_tran.m` | `A=sin_tran(dim)` | Single transition operators, spanning the space of matri- ces of the specified dimension. The set is returned as a cell |
| `kernel/operators/stevens.m` | `S=stevens(mult,k,q)` | Extended Stevens operators. Syntax: S=stevens(mult,k,q) Parameters: mult -multiplicity of the spin in question k -Steven |
| `kernel/operators/superop.m` | `A=superop(spin_system,opspec,side)` | Sided product superoperator in the spherical tensor basis set. Returns superoperators corresponding to right or left mul |
| `kernel/operators/twospinist.m` | `T=twospinist(spin_system,spin_a,spin_b,indices,type)` | Two-spin irreducible spherical tensor operators. Syntax: T=twospinist(spin_system,spin_a,spin_b,indices,type) Parameters |
| `kernel/operators/unit_oper.m` | `A=unit_oper(spin_system)` | Returns a unit operator in the current formalism and basis. The operator has dimension equal to the basis size in sphten |
| `kernel/operators/weyl.m` | `A=weyl(nlevels)` | Weyl boson operators (sparse, see below for normalisa- tion convention) for a bosonic mode with a user-speci- fied popul |
| `kernel/optimcon/alpha_conds.m` | `test=alpha_conds(test_type,alpha,fx_0,fx_1,gfx_0,gfx_1,dir,spin_system)` | Applies one of the line search acceptance tests used by the brac- keting and sectioning routines in constrained optimisa |
| `kernel/optimcon/aux_mat.m` | `[auxm_l,auxm_r]=aux_mat(drifts,controls,cc_comm_idx,cc_comm,dt,cL,cR,k,j)` | Builds auxiliary matrices for the calculation of the directional derivatives of the trapezium product quadrature propaga |
| `kernel/optimcon/bfgs.m` | `H=bfgs(dx_hist,dg_hist,g)` | Calculates a BFGS approximation to the Newton-Raphson search direction for maximising a function using past gradients to |
| `kernel/optimcon/bfgs_upd.m` | `H=bfgs_upd(H,dx,dg)` | Performs a one-step BFGS Hessian update for maximisation using the argument and gradient increments from the previous st |
| `kernel/optimcon/bracketing.m` | `[a,b,alpha,fx,gfx,next_act,data]=bracketing(cost_function,alpha,dir,x_0,fx_0,gfx_0,data,spin_system)` | Expands a trial step into a bracket that contains an acceptable line search point, or accepts the step directly if the W |
| `kernel/optimcon/bss_ops.m` | `resp_ops=bss_ops(spin_system,channels,carrier_frq)` | Bloch-Siegert response operators for the optimal control module. For each control channel, returns the operator whose co |
| `kernel/optimcon/ctrl_trajan.m` | `ctrl_trajan(spin_system,waveform,traj_data,fidelities)` | Diagnostic plotting function for optimal control module. Plots trajectory and control pulse analysis. Syntax: ctrl_traja |
| `kernel/optimcon/cubic_interp.m` | `[alpha,fx]=cubic_interp(end_a,end_b,alpha_a,alpha_b,f_a,dir_der_a,f_b,dir_der_b)` | Finds the extremum of a cubic interpolant built from function values and directional derivatives at two points and retur |
| `kernel/optimcon/distortions/amp_root.m` | `[w,J]=amp_root(w,sat_lvls,s)` | Amplifier compression distortion model. Applies a saturating root-sigmoidal distortion to the waveform amplitude: y=x/(1 |
| `kernel/optimcon/distortions/amp_tanh.m` | `[w,J]=amp_tanh(w,sat_lvls)` | Amplifier compression distortion model. Applies a saturating hyperbolic tangent distortion to the user-supplied waveform |
| `kernel/optimcon/distortions/firf.m` | `[w,J]=firf(w,ker)` | Applies an FIR convolution filter to a Spinach optimal control module waveform. Treats odd rows of multi-row waveform ar |
| `kernel/optimcon/distortions/kernelest.m` | `h=kernelest(x,y,ker_len,method,align,lambda)` | FIR convolution kernel estimation from input and output signal samples. Syntax: h=kernelest(x,y,ker_len,method,align,lam |
| `kernel/optimcon/distortions/no_dist.m` | `[w,J]=no_dist(w)` | A distortion function that applies no distortion and therefore has a unit Jacobian. Syntax: [w,J]=no_dist(w) Parameters: |
| `kernel/optimcon/distortions/non_orth.m` | `[w,J]=non_orth(w,xy_ang)` | Non-orthogonal channel distortion model. Treats odd rows of multi-row waveform arrays as in-phase channels, and even row |
| `kernel/optimcon/distortions/spf.m` | `[w,J]=spf(w,p)` | Applies a discrete single-pole filter: Y(n)=(1-p)*X(n)+p*Y(n-1) to a Spinach optimal control module waveform. Treats odd |
| `kernel/optimcon/distortions/szf.m` | `[w,J]=szf(w,z)` | Applies a discrete single-zero filter: Y(k)=X(k)/(1-z)-z*X(k-1)/(1-z); to a Spinach optimal control module waveform. Tre |
| `kernel/optimcon/drifts.m` | `[drifts,spc_dim]=drifts(spin_system,context,parameters,assumptions)` | Returns a cell array of drift Liouvillians suitable for the control.drifts variable in ensemble control optimisations. S |
| `kernel/optimcon/ens_block.m` | `[traj,fid,grad,hess]=ens_block(spin_system,drifts,control,block,waveform,n_outputs)` | Fidelity, gradient, and Hessian contributions of one block of ensemble cases, evaluated on the parallel pool worker that h |
| `kernel/optimcon/ens_catalog.m` | `[catalog,ens_sizes]=ens_catalog(control)` | Ensemble case catalog for optimal control problems. Enumerates the Cartesian product of the state-target pairs, the drif |
| `kernel/optimcon/ensemble.m` | `[traj_data,fidelity,gradient,hessian]=ensemble(waveform,spin_system)` | A parallel wrapper around GRAPE that enables ensemble optimal control optimisations. This function handles systems with |
| `kernel/optimcon/fapt2sfo.m` | `[wave,dt,time_grid]=fapt2sfo(fapt,time_grid)` | Converts a freq-ampl-phase-time specification of a pulse sequ- uence into the corresponding single frequency origin wave |
| `kernel/optimcon/fmaxnewton.m` | `[x,data]=fmaxnewton(spin_system,cost_function,guess)` | Finds a local maximum of a function of several variables using Newton and quasi-Newton algorithms. Syntax: [x,data]=fmax |
| `kernel/optimcon/grape_hilb.m` | `[traj_data,fidelity,grad,hess]=grape_hilb(spin_system,drifts,controls,waveform,rho_init,rho_targ,fidelity_type)` | Gradient Ascent Pulse Engineering (GRAPE) objective function, gradient and Hessian. Propagates the system through a user |
| `kernel/optimcon/grape_liouv.m` | `[traj_data,fidelity,grad,hess]=grape_liouv(spin_system,drifts,controls,waveform,rho_init,rho_targ,fidelity_type)` | Gradient Ascent Pulse Engineering (GRAPE) objective function, gradient and Hessian. Propagates the system through a user |
| `kernel/optimcon/hess_reorder.m` | `hess=hess_reorder(hess,K,N)` | The waveforms on different channels are assumed to be stored in the rows of the input array. The Hessian elements corres |
| `kernel/optimcon/hessreg.m` | `[H,data]=hessreg(spin_system,H,g,data)` | RFO regularisation for Newton-Raphson Hessian and gradient pairs. Syntax: [H,data]=hessreg(spin_system,H,g,data) Paramet |
| `kernel/optimcon/inst_freq.m` | `freq=inst_freq(signal,dt,npoints,poly_order,amp_tol)` | Instantaneous frequency trajectory from a complex time-domain signal by regularised phase differentiation. Syntax: freq= |
| `kernel/optimcon/lbfgs.m` | `direction=lbfgs(dx_hist,dg_hist,g)` | Calculates an approximation to the Newton-Raphson search direction for maximising a function using past gradients to bui |
| `kernel/optimcon/mas_drifts.m` | `drifts=mas_drifts(spin_system,parameters)` | Rotor-phase-resolved MAS drift Hamiltonians over powder orientations and initial rotor phases. |
| `kernel/optimcon/objeval.m` | `[data,fx,grad,hess]=objeval(x,objfun_handle,data,spin_system)` | Calls and collect the correct amount of outputs from an objective function -used by optimisation routines. Syntax: [data |
| `kernel/optimcon/optimcon.m` | `spin_system=optimcon(spin_system,control)` | Validates optimal control options and updates the spin system object. Syntax: spin_system=optimcon(spin_system,control) |
| `kernel/optimcon/penalty.m` | `[pen_term,pen_grad,pen_hess]=penalty(wf,type,fb,cb)` | Penalty terms for the Optimal Control module. Returns the penalty function and its gradient for the waveform, which shou |
| `kernel/optimcon/sectioning.m` | `[alpha,fx_1,gfx_1,exitflag,data]=sectioning(cost_function,a,b,x_0,fx_0,gfx_0,dir,data,spin_system)` | Refines a previously found step bracket by repeated cubic interpolation until a step satisfying Wolfe tests is found or |
| `kernel/optimcon/tgrape.m` | `[fidelity,grad]=tgrape(spin_system,drift,controls,waveform,dt_grid,time_unit,rho_init,rho_targ)` | A special case of Gradient Ascent Pulse Engineering (GRAPE) objective function and gradient with respect to the vector o |
| `kernel/optimcon/trapdiff.m` | `[DL,DR]=trapdiff(spin_system,Hd,Hc,dt,cL,cR)` | Directional derivatives for the trapezium product quadrature publi- shed by Iserles and Norsett (see Corollary 3.3) in T |
| `kernel/optimcon/wrappers/grape_coop.m` | `[traj_data,fidelity,gradient]=grape_coop(phi_profile,spin_system)` | Pairs of cooperative pulses that may be used as components of a phase cycle. The pulses are designed to produce as much |
| `kernel/optimcon/wrappers/grape_curv.m` | `[traj_data,fidelity,df_du]=grape_curv(waveform_u,u2x,dx_du,spin_system)` | Cost function for optimal control using the GRAPE algorithm. Returns fidelity and gradient for a given waveform, specifi |
| `kernel/optimcon/wrappers/grape_phase.m` | `[traj_data,fidelity,gradient,hessian]=grape_phase(phi_profile,spin_system)` | Cost function for optimal control using the GRAPE algorithm. Returns fidelity, gradient and Hessian for a given waveform |
| `kernel/optimcon/wrappers/grape_xy.m` | `[traj_data,fidelity,grad,hess]=grape_xy(waveform,spin_system)` | Cost function for optimal control using the GRAPE algorithm. Returns fidelity, gradient and hessian for a given waveform |
| `kernel/orientation.m` | `H=orientation(Q,euler_angles)` | Anisotropic part of the Hamiltonian for a specific spin system orientation. Syntax: H=orientation(Q,euler_angles) Argume |
| `kernel/overloads/@cell/blkdiag.m` | `C=blkdiag(A,B)` | Block-diagonal cell array from two cell arrays, all other elements are set to empty cells. Syntax: C=blkdiag(A,B) Parame |
| `kernel/overloads/@cell/complex.m` | `A=complex(A)` | A shorthand for making all elements of a cell array complex. Syntax: A=complex(A) Parameters: A - a cell array of numeri |
| `kernel/overloads/@cell/inflate.m` | `A=inflate(A)` | A shorthand for inflating cell arrays of polyadics; this inflates every component of the cell array. Syntax: A=inflate(A |
| `kernel/overloads/@cell/minus.m` | `C=minus(A,B)` | Subtracts cell arrays element-by-element. Syntax: A=minus(A,B) Parameters: A,B -cell arrays of identical topology Output |
| `kernel/overloads/@cell/mtimes.m` | `C=mtimes(A,B)` | Multiplies all entries of a cell array by a user-specified matrix. Syntax: C=mtimes(A,B) Parameters: A -a matrix or a ce |
| `kernel/overloads/@cell/plus.m` | `C=plus(A,B)` | Adds cell arrays element-by-element. Syntax: A=plus(A,B) Parameters: A,B -cell arrays of identical topology Outputs: C - |
| `kernel/overloads/@cell/times.m` | `C=times(A,B)` | Multiplies all entries of a cell array by a user-specified scalar or a matching dimension numeric array. Syntax: C=times |
| `kernel/overloads/@cell/totsum.m` | `S=totsum(A)` | A sum across all dimensions of a cell array. Syntax: S=totsum(A) Parameters: A -a cell array of numerical objects Output |
| `kernel/overloads/@double/inflate.m` | `A=inflate(A)` | A dummy function that mirrors the inflate() command for polyadics, but is not supposed to do anything to proper numerica |
| `kernel/overloads/@opium/full.m` | `M=full(M)` | Converts an OPIUM object into the full scaled unit matrix that it represents. Syntax: M=full(M) Parameters: M -an OPIUM |
| `kernel/overloads/@opium/kron.m` | `c=kron(a,b)` | Kronecker products involving an OPIUM object. Syntax: c=kron(a,b) Parameters: a,b -Kronecker operands, can be matrices o |
| `kernel/overloads/@opium/mtimes.m` | `c=mtimes(a,b)` | Matrix products involving an OPIUM object. Syntax: c=mtimes(a,b) Parameters: a,b -opia or numerical arrays Outputs: c -m |
| `kernel/overloads/@opium/opium.m` | `M=opium(dim,coeff)` | Object Pretending It is a Unit Matrix (OPIUM). Syntax: M=opium(dim,coeff) Parameters: dim -dimension of the unit matrix |
| `kernel/overloads/@opium/size.m` | `varargout=size(op,dim)` | The size of the matrix represented by the OPIUM. Syntax: answer=size(op,dim) Parameters: op -an opium object dim -option |
| `kernel/overloads/@polyadic/allfinite.m` | `answ=allfinite(p)` | Returns true if none of the elements of the polyadic are Inf or NaN. Syntax: answ=allfinite(p) Parameters: p -a polyadic |
| `kernel/overloads/@polyadic/ctranspose.m` | `p=ctranspose(p)` | Computes the Hermitian conjugate of a matrix in a polyadic representation. Syntax: p=ctranspose(p) |
| `kernel/overloads/@polyadic/full.m` | `answer=full(p)` | Converts a polyadic representation of a matrix into a full mat- rix. Syntax: answer=full(p) The function opens up all th |
| `kernel/overloads/@polyadic/gpuArray.m` | `p=gpuArray(p)` | Uploads all components of a polyadic object to the GPU. The object still looks like a polyadic to Matlab, but all of its |
| `kernel/overloads/@polyadic/inflate.m` | `answer=inflate(p)` | Converts a polyadic representation of a matrix into a sparse mat- rix. Syntax: answer=inflate(p) The function opens up a |
| `kernel/overloads/@polyadic/isempty.m` | `answer=isempty(p)` | Returns true for polyadics that represent a matrix with a zero dimension. Syntax: answer=isempty(p) Parameters: p -a pol |
| `kernel/overloads/@polyadic/isreal.m` | `answer=isreal(p)` | Returns true if the polyadic representation is real. Syntax: answer=isreal(p) Parameters: p -a polyadic object Outputs: |
| `kernel/overloads/@polyadic/kron.m` | `c=kron(a,b)` | Kronecker product function for polyadics. Syntax: c=kron(a,b) Parameters: a,b -polyadic or numeric objects Outputs: c -p |
| `kernel/overloads/@polyadic/minus.m` | `a=minus(a,b)` | Polyadic subtraction operation. Does not perform the actual sub- traction, but instead stores the operands as a sum of u |
| `kernel/overloads/@polyadic/mtimes.m` | `C=mtimes(A,B)` | Performs multiplications involving polyadics. Syntax: C=mtimes(A,B) Parameters: A,B -a polyadic or a numerical array Out |
| `kernel/overloads/@polyadic/nnz.m` | `answer=nnz(p)` | Number of non-zeroes in all kernels of the polyadic. Syntax: answer=nnz(p) Parameters: p -a polyadic object Outputs: ans |
| `kernel/overloads/@polyadic/plus.m` | `c=plus(a,b)` | Polyadic addition operation. Does not perform the actual additi- on, but instead stores the operands as a sum of unopene |
| `kernel/overloads/@polyadic/polyadic.m` | `p=polyadic(cores)` | Creates an object of a polyadic class. Syntax: p=polyadic(cores) A polyadic is a matrix formed by a Kronecker product, w |
| `kernel/overloads/@polyadic/prefix.m` | `p=prefix(a,p)` | Adds prefix matrices to a polyadic. Anything the polyadic multiplies will subsequently be multiplied by the prefix matri |
| `kernel/overloads/@polyadic/simplify.m` | `p=simplify(p)` | Simplifies the structure of the polyadic object by reordering buffers, dropping inconsequential terms, and flattening ne |
| `kernel/overloads/@polyadic/size.m` | `varargout=size(p,dim)` | Returns the size of the matrix represented by the polyadic. Syntax: answer=size(p,dim) Parameters: p -a polyadic object |
| `kernel/overloads/@polyadic/suffix.m` | `p=suffix(p,a)` | Adds suffix matrices to a polyadic. Anything the polyadic multiplies will first be multiplied by the suffix matri- ces. |
| `kernel/overloads/@polyadic/transpose.m` | `p=transpose(p)` | Computes the transpose of a matrix in a polyadic representa- tion. Syntax: p=transpose(p) |
| `kernel/overloads/@polyadic/validate.m` | `validate(p)` | Checks the internal structure of a polyadic object and throws an error if the object does not meet expectations. Syntax: |
| `kernel/overloads/@rcv/ctranspose.m` | `A=ctranspose(A)` | Returns the conjugate transpose of an RCV sparse matrix. Syntax: A=ctranspose(A) Parameters: A -an RCV sparse matrix Out |
| `kernel/overloads/@rcv/full.m` | `A=full(A)` | Converts an RCV sparse matrix into a full matrix. Syntax: A=full(A) Parameters: A -an RCV sparse matrix Outputs: A -a fu |
| `kernel/overloads/@rcv/gather.m` | `A=gather(A)` | Gathers an RCV sparse matrix from GPU. Syntax: A=gather(A) Parameters: A -an RCV sparse matrix Outputs: A -the same matr |
| `kernel/overloads/@rcv/gpuArray.m` | `obj=gpuArray(obj)` | Transfers an RCV sparse matrix to the GPU. Syntax: obj=gpuArray(obj) Parameters: obj -an RCV sparse matrix Outputs: obj |
| `kernel/overloads/@rcv/horzcat.m` | `A=horzcat(varargin)` | Horizontal concatenation for RCV sparse matrices. Syntax: A=horzcat(A,B) Parameters: A -left RCV sparse matrix B -right |
| `kernel/overloads/@rcv/kron.m` | `C=kron(A,B)` | Kronecker product between two RCV sparse matrices. Syntax: C=kron(A,B) Parameters: A -left RCV sparse matrix B -right RC |
| `kernel/overloads/@rcv/minus.m` | `A=minus(A,B)` | Subtracts one RCV object from another. Syntax: A=minus(A,B) Parameters: A -left operand B -right operand Outputs: A -res |
| `kernel/overloads/@rcv/mtimes.m` | `C=mtimes(A,B)` | Multiplication for RCV sparse matrices. Syntax: C=mtimes(A,B) Parameters: A -left operand B -right operand Outputs: C -p |
| `kernel/overloads/@rcv/plus.m` | `C=plus(A,B)` | Adds things to RCV sparse matrices. Syntax: C=plus(A,B) Parameters: A -left operand B -right operand Outputs: C -sum A+B |
| `kernel/overloads/@rcv/rcv.m` | `obj=rcv(varargin)` | Creates an RCV (row-column-value storage) sparse matrix. Syntax: obj=rcv(M) obj=rcv(dim1,dim2) obj=rcv(R,C,V,dim1,dim2) |
| `kernel/overloads/@rcv/rdivide.m` | `A=rdivide(A,k)` | Divides an RCV sparse matrix by a numeric scalar. Syntax: A=rdivide(A,k) Parameters: A -RCV sparse matrix k -numeric sca |
| `kernel/overloads/@rcv/size.m` | `[s,ncols]=size(A,dim)` | Returns the size of an RCV sparse matrix. Syntax: s=size(A,dim) [s,ncols]=size(A) Parameters: A -RCV sparse matrix dim -optional dimension |
| `kernel/overloads/@rcv/sparse.m` | `A=sparse(A)` | Converts an RCV sparse matrix into a Matlab sparse matrix. Syntax: A=sparse(A) Parameters: A -RCV sparse matrix Outputs: |
| `kernel/overloads/@rcv/spy.m` | `spy(A)` | Plots the sparsity pattern of an RCV matrix. Syntax: spy(A) Parameters: A -RCV sparse matrix Outputs: produces a sparsit |
| `kernel/overloads/@rcv/times.m` | `C=times(A,B)` | Multiplies an RCV sparse matrix by a numeric scalar, in either operand order. Syntax: C=times(A,B) Parameters: A,B -an R |
| `kernel/overloads/@rcv/transpose.m` | `A=transpose(A)` | The transpose of an RCV sparse matrix. Syntax: A=transpose(A) Parameters: A -an RCV sparse matrix Outputs: A -transposed |
| `kernel/overloads/@rcv/vertcat.m` | `A=vertcat(varargin)` | Vertical concatenation for RCV sparse matrices. Syntax: A=vertcat(A,B) Parameters: A -top RCV matrix B -bottom RCV matri |
| `kernel/overloads/@struct/mtimes.m` | `str_out=mtimes(M,str_in)` | Multiplies all entries of a structure by a user-specified mat- rix. Nested structures are processed recursively. Syntax: |
| `kernel/overloads/@struct/plus.m` | `str3=plus(str1,str2)` | Adds corresponding fields of two structures. Nested structu- res are processed recursively. Syntax: str3=plus(str1,str2) |
| `kernel/overloads/@ttclass/amensolve.m` | `x=amensolve(A,y,tol,opts,x0) % #NORMOK` | Solves the linear system Ax=y using the AMEn iteration. Syntax: x=amensolve(A,y,tol,opts,x0) Parameters: A -ttclass repr |
| `kernel/overloads/@ttclass/amensum.m` | `y=amensum(x,tol,opts)` | Sums buffered tensor trains in a single tensor train using AMEn algorithm. Syntax: y=amensum(x,tol,opts) Parameters: x - |
| `kernel/overloads/@ttclass/clearcoeff.m` | `tt=clearcoeff(tt)` | Absorbs physical coefficients into tensor train cores without changing the value represented by the tensor train. Syntax |
| `kernel/overloads/@ttclass/conj.m` | `tt=conj(tt)` | Conjugates all core elements and coefficients of a tensor train object. Syntax: tt=conj(tt) Parameters: tt -tensor train |
| `kernel/overloads/@ttclass/ctranspose.m` | `ttrain=ctranspose(ttrain)` | Computes a Hermitian conjugate of a matrix in a tensor train representation. Syntax: ttrain=ctranspose(ttrain) Parameter |
| `kernel/overloads/@ttclass/diag.m` | `tt=diag(tt)` | Mimics the diag behaviour for tensor train matrix. Syntax: tt=diag(tt) Parameters: tt -a tensor train representation of |
| `kernel/overloads/@ttclass/dot.m` | `c=dot(a,b)` | Dot product of TT representations of matrices. Syntax: c=dot(a,b) Parameters: a,b -tensor train objects representing num |
| `kernel/overloads/@ttclass/full.m` | `answer=full(ttrain)` | Converts a tensor train representation of a matrix into a matrix. Syntax: answer=full(ttrain) Parameters: ttrain -tensor |
| `kernel/overloads/@ttclass/hdot.m` | `c=hdot(a,b)` | Hadamard dot product between two tensor train matrices. Syntax: c=hdot(a,b) Parameters: a,b -tensor train objects repres |
| `kernel/overloads/@ttclass/ismatrix.m` | `answer=ismatrix(tt)` | Returns TRUE for non-empty tensor train objects. Syntax: answer=ismatrix(tt) Parameters: tt -tensor train object Outputs |
| `kernel/overloads/@ttclass/isnumeric.m` | `answer=isnumeric(tt)` | Returns TRUE for non-empty tensor train objects. Syntax: answer=isnumeric(tt) Parameters: tt -tensor train object Output |
| `kernel/overloads/@ttclass/isreal.m` | `answer=isreal(tt)` | Returns TRUE for real-valued tensor train objects. Syntax: answer=isreal(tt) Parameters: tt -tensor train object Outputs |
| `kernel/overloads/@ttclass/kron.m` | `c=kron(a,b)` | Kronecker product of two matrices in a tensor train format. Syntax: c=kron(a,b) Parameters: a,b -tensor train objects Ou |
| `kernel/overloads/@ttclass/mean.m` | `answer=mean(ttrain,dim)` | Mean of elements of a tensor train representation of a matrix. Syntax: answer=mean(ttrain,dim) Parameters: ttrain -a ten |
| `kernel/overloads/@ttclass/minus.m` | `a=minus(a,b)` | Tensor train subtraction operation. Does not perform the actual subtraction but instead concatenates the operands until |
| `kernel/overloads/@ttclass/mldivide.m` | `x=mldivide(A,y)` | Solves a linear system with tensor train objects. Syntax: x=mldivide(A,y) Parameters: A -ttclass matrix y -ttclass vecto |
| `kernel/overloads/@ttclass/mrdivide.m` | `a=mrdivide(a,b)` | Divides a tensor train object by a scalar. Syntax: c=mrdivide(a,b) Parameters: a -tensor train object b -a scalar Output |
| `kernel/overloads/@ttclass/mtimes.m` | `c=mtimes(a,b)` | Performs tensor train multiplication followed by a shrink. Syntax: c=mtimes(a,b) Parameters: a -a scalar or a tensor tra |
| `kernel/overloads/@ttclass/nnz.m` | `answer=nnz(ttrain)` | Counts non-zero elements in all cores of a tensor train. Syntax: answer=nnz(ttrain) Parameters: ttrain -tensor train obj |
| `kernel/overloads/@ttclass/norm.m` | `ttnorm=norm(ttrain,norm_type) %#NORMOK` | Computes the norm of the matrix represented by a tensor train. Syntax: ttnorm=norm(ttrain,norm_type) Parameters: ttrain |
| `kernel/overloads/@ttclass/numel.m` | `n=numel(tt)` | Number of elements in the matrix represented by a tensor train. Syntax: n=numel(tt) Parameters: tt -tensor train object |
| `kernel/overloads/@ttclass/pack.m` | `ttout=pack(tt)` | This subroutine packs all trains from the addition buffer into a single tensor train, but does not perform the recom- pr |
| `kernel/overloads/@ttclass/plus.m` | `a=plus(a,b)` | Tensor train addition operation. Does not perform the actual addition, but instead concatenates the operands until such |
| `kernel/overloads/@ttclass/rand.m` | `tt=rand(tt,ttrank)` | Generates a tensor train representation of a matrix with random tensor train cores, same physical index topology as the |
| `kernel/overloads/@ttclass/ranks.m` | `ttranks=ranks(ttrain)` | Returns the bond dimensions of a tensor train. Syntax: ttranks=ranks(ttrain) Parameters: ttrain -a tensor train object O |
| `kernel/overloads/@ttclass/rdivide.m` | `a=rdivide(a,b)` | Divides a tensor train object by a scalar. Syntax: c=rdivide(a,b) Parameters: a -a ttclass object b -a numeric scalar Ou |
| `kernel/overloads/@ttclass/revert.m` | `tt=revert(tt)` | Applies a bit-revert permutation to a tensor train operator by reversing the core order and swapping bond indices. Synta |
| `kernel/overloads/@ttclass/shrink.m` | `ttrain=shrink(ttrain)` | Approximates a given tensor train with lower TT-ranks. Syntax: ttrain=shrink(ttrain) Parameters: ttrain -a tensor train |
| `kernel/overloads/@ttclass/size.m` | `varargout=size(tt,dim)` | Returns the size of the matrix represented by a tensor train. The output mimics Matlab's size function. Syntax: [m,n]=si |
| `kernel/overloads/@ttclass/sizes.m` | `modesizes=sizes(tt)` | Returns mode sizes (physical dimensions of each core) of a tensor train. Syntax: modesizes=sizes(tt) Parameters: tt -ten |
| `kernel/overloads/@ttclass/subsref.m` | `answer=subsref(ttrain,reference)` | Dot and bracket property specifications for the tensor train class. Syntax: answer=subsref(ttrain,reference) Parameters: |
| `kernel/overloads/@ttclass/sum.m` | `answer=sum(ttrain,dim)` | Sum of elements of a tensor train representation of a matrix. Syntax: answer=sum(ttrain,dim) Parameters: ttrain -tensor |
| `kernel/overloads/@ttclass/trace.m` | `tttrace=trace(tt)` | Computes the trace of a tensor train operator. Syntax: tttrace=trace(tt) Parameters: tt -tensor train operator Outputs: |
| `kernel/overloads/@ttclass/transpose.m` | `ttrain=transpose(ttrain)` | Transposes a tensor without complex conjugation. Syntax: ttrain=transpose(ttrain) Parameters: ttrain -tensor train repre |
| `kernel/overloads/@ttclass/truncate.m` | `ttout=truncate(tt)` | Performs right-to-left SVD recompression for a tensor train. This should not be called directly, use shrink.m instead. S |
| `kernel/overloads/@ttclass/ttclass.m` | `tt=ttclass(coeff,kronterms,tolerance)` | Creates an object of a tensor train class. A tensor train is a type of un-opened Kronecker product that behaves as a mat |
| `kernel/overloads/@ttclass/ttort.m` | `[tt,lognrm]=ttort(tt,direct)` | Performs TT-orthogonalisation for a tensor train (or for each tensor train in a buffered sum). Syntax: [tt,lognrm]=ttort |
| `kernel/overloads/@ttclass/unit_like.m` | `A=unit_like(A)` | Returns a unit object of the same type as whatever is supplied. Syntax: A=unit_like(A) Parameters: A -a full or sparse s |
| `kernel/overloads/@ttclass/vec.m` | `A=vec(A)` | Stretches arrays into vectors -useful for situations when the stand- ard (:) syntax is not available. Syntax: A=vec(A) P |
| `kernel/overloads/save_anyway.m` | `save_anyway(file_name,variable)` | A wrapper intended to trick SPMD blocks into saving data. Can only save one variable at a time, its name in the mat file |
| `kernel/plotting/axis_1d.m` | `[ax,ax_label]=axis_1d(spin_system,parameters)` | Generates axis ticks for plotting 1D spectra. Syntax: [ax,ax_label]=axis_1d(spin_system,parameters) Parameters: paramete |
| `kernel/plotting/bloch_axis.m` | `[ax,ay,az]=bloch_axis(x,y,z)` | Reconstructs the instantaneous Bloch equation rotation axis of from a 3D magnetisation trajectory. Syntax: [ax,ay,az]=bl |
| `kernel/plotting/bwr_cmap.m` | `cmap=bwr_cmap()` | Blue -> White -> Red colour map with 255 points and white colour corresponding to zero. Syntax: cmap=bwr_cmap() The outp |
| `kernel/plotting/contspacing.m` | `[all_conts,pos_conts,neg_conts]=contspacing(smax,smin,delta,k,signs,ncont)` | Non-linear adaptive contour spacing. Useful for NMR data where small cross-peaks must be adequately contoured next to la |
| `kernel/plotting/crop_2d.m` | `[spec,parameters]=crop_2d(spin_system,spec,parameters,crop_ranges)` | Crops 2D spectra to user-specified ranges (in ppm), respecting the digital resolution. Syntax: [spec,parameters]=crop_2d |
| `kernel/plotting/cst_display.m` | `cst_display(props,atoms,scaling,conmatrix,options)` | Draws shielding tensors and their eigensystems. Two styles are implemented: A. Ellipsoids (symmetric tensors only): 1. A |
| `kernel/plotting/cylgrid.m` | `cylgrid(zmin,zmax,rmax)` | Draws a cylindrical grid with 10% spacing added around the indicated data extent values. Syntax: cylgrid(zmin,zmax,rmax) |
| `kernel/plotting/efg_display.m` | `efg_display(props,atoms,scaling,conmatrix,options)` | Electric field gradient tensors and their eigensystems. Two styles are implemented: A. Ellipsoids (symmetric tensors onl |
| `kernel/plotting/fft_freq_axis.m` | `[f_shift,f,df]=fft_freq_axis(npts,dt,zf)` | Frequency axis for FFT with optional zero-filling. Syntax: [f_shift,f,df,nfft]=fft_freq_axis(npts,dt,zf) Parameters: npt |
| `kernel/plotting/fig2tiles.m` | `[fig_obj,tile_obj]=fig2tiles(fig_files,fig_size)` | Combines Matlab figure files into a single tiled figure. Syntax: [fig_obj,tile_obj]=fig2tiles(fig_files,fig_size) Parame |
| `kernel/plotting/ft_axis.m` | `ax=ft_axis(offset,sweep,npoints)` | Fourier transform axis ticks generator that accounts for the periodicity and correctly folds the edge frequency. Syntax: |
| `kernel/plotting/hfc_display.m` | `hfc_display(props,atoms,scaling,conmatrix,options)` | Draws hyperfine tensors and their eigensystems. Two styles are implemented: A. Ellipsoids (symmetric tensors only): 1. A |
| `kernel/plotting/ifft_time_axis.m` | `[t_shift,t,dt]=ifft_time_axis(npts,df,zf)` | Time axis for IFFT with optional zero-filling. Syntax: [t_shift,t,dt,nifft]=ifft_time_axis(npts,df,zf) Parameters: npts |
| `kernel/plotting/int_2d.m` | `int_2d(spin_system,spectrum,parameters,ncont,delta,k,ncol,m,signs,filename)` | Contour plotting utility with non-linear adaptive contour spacing and 2D integration using mouse or an interval file. Sy |
| `kernel/plotting/kbox.m` | `kbox() % #NGRUM` | Creates a tickless boxed frame around the current axes using ordinary line objects that live in the same data space as t |
| `kernel/plotting/kcolourbar.m` | `kcolourbar(x)` | House style settings for Matlab figures; a product of much experience with academic publication aesthetics. Syntax: kcol |
| `kernel/plotting/kfigure.m` | `handle=kfigure(varargin)` | Resets the stupid ass figure defaults in R2025a and later back to sensible values. |
| `kernel/plotting/kgrid.m` | `kgrid()` | A replacement for the 'grid' command in Matlab that produces grey (rather than black-and-transparent) grid lines that ar |
| `kernel/plotting/klegend.m` | `leg_obj=klegend(varargin)` | House style settings for Matlab figures; a product of much experience with academic publication aesthetics. Syntax: leg_ |
| `kernel/plotting/kletter.m` | `kletter(letter_label)` | Draws an academic journal style letter label in the top left corner of the current axis set. The label is placed inside the |
| `kernel/plotting/ksgtitle.m` | `ksgtitle(x)` | House style settings for Matlab figures; a product of much experience with academic publication aesthetics. Syntax: ksgt |
| `kernel/plotting/ktitle.m` | `ktitle(x)` | House style settings for Matlab figures; a product of much experience with academic publication aesthetics. Syntax: ktit |
| `kernel/plotting/kxlabel.m` | `kxlabel(varargin)` | House style settings for Matlab figures; a product of much experience with academic publication aesthetics. Syntax: kxla |
| `kernel/plotting/kxtickfix.m` | `kxtickfix()` | Switches X axis tick labels to engineering notation by setting the numeric-ruler exponent to a multiple of three and ins |
| `kernel/plotting/kylabel.m` | `kylabel(varargin)` | House style settings for Matlab figures; a product of much experience with academic publication aesthetics. Syntax: kyla |
| `kernel/plotting/kytickfix.m` | `kytickfix()` | Switches Y axis tick labels to engineering notation by setting the numeric-ruler exponent to a multiple of three and ins |
| `kernel/plotting/kzlabel.m` | `kzlabel(varargin)` | House style settings for Matlab figures; a product of much experience with academic publication aesthetics. Syntax: kzla |
| `kernel/plotting/kztickfix.m` | `kztickfix()` | Switches Z axis tick labels to engineering notation by setting the numeric-ruler exponent to a multiple of three and ins |
| `kernel/plotting/molplot.m` | `molplot(xyz,conmatrix)` | Plots a stick representation of a molecule from Cartesian coordinates supplied. Syntax: molplot(xyz,conmatrix) Parameter |
| `kernel/plotting/mri_2d_plot.m` | `mri_2d_plot(mri,parameters,method)` | 2D MRI image plotting. Syntax: mri_2d_plot(mri,parameters,method) Parameters: method -'image' uses gradient information |
| `kernel/plotting/plot_1d.m` | `plot_1d(spin_system,spectrum,parameters,varargin)` | 1D plotting utility. Syntax: plot_1d(spin_system,spectrum,parameters,varargin) Parameters: spectrum a column vector cont |
| `kernel/plotting/plot_2d.m` | `[axis_f1,axis_f2,spectrum]=plot_2d(spin_system,spectrum,parameters,ncont,delta,k,ncol,m,signs)` | Contour plotting utility with non-linear adaptive contour spacing. The function is useful for NMR data where small cross |
| `kernel/plotting/plot_3d.m` | `plot_3d(spin_system,spectrum,parameters,nsurf,delta,k,signs)` | Volume isosurface plotting utility with non-linear adaptive surface spacing. The function plots the 3D volume and the th |
| `kernel/plotting/plot_uf.m` | `plot_uf(spin_system,spectrum_uf,parameters)` | Plotting utility for ultrafast constant-time 2D pulse sequences. Syntax: plot_uf(spin_system,spectrum_uf,parameters) Par |
| `kernel/plotting/scale_figure.m` | `scale_figure(by)` | Scales the current figure from the default size by the factors provided by the user. Syntax: scale_figure(by) Parameters |
| `kernel/plotting/slice_2d.m` | `slice_2d(spin_system,spectrum,parameters,ncont,delta,k,ncol,m,signs)` | Contour plotting utility with non-linear adaptive contour spacing and 1D slice extraction using mouse. Syntax: slice_2d( |
| `kernel/plotting/stack_2d.m` | `stack_2d(spin_system,spectrum,parameters,stack_dim,alpha_fun)` | Stack plotting utility for 2D NMR spectra. Syntax: stack_2d(spin_system,spectrum,parameters,stack_dim,alpha_fun) Paramet |
| `kernel/plotting/sweep2ticks.m` | `axis_hz=sweep2ticks(offs,sweep,npoints)` | Converts offset-sweep-npoints specification into axis ticks in Hz. The function returns the frequency axis of the spectr |
| `kernel/plotting/volplot.m` | `volplot(data_cube,axis_ranges,clip_ranges)` | Volumetric 3D plot function for scalar fields. Sign is mapped into colour and amplitude into opacity. Separate scaling f |
| `kernel/plotting/write_movie.m` | `write_movie(file_name)` | Orbits the camera around a 3D plot and writes a correctly looping 359-frame movie. Intended for generating presenta- tio |
| `kernel/plotting/zoom_3d.m` | `[density,ext]=zoom_3d(density,ext,zoom_ranges)` | Zooms a 3D data cube to the fractional limits specified by the user. Syntax: [density,ext]=zoom_3d(density,ext,zoom_rang |
| `kernel/ppower.m` | `P=ppower(spin_system,P,N)` | Computes integer propagator powers via an efficient powers-of-two based strategy. Syntax: P=ppower(spin_system,P,N) Para |
| `kernel/propagator.m` | `P=propagator(spin_system,L,timestep)` | Calculates exponential propagator exp(-1i*L*timestep) using scaled and squared Taylor series method. Syntax: P=propagato |
| `kernel/pulses/bloch_siegert.m` | `[ctrl_opers,ctrl_coefs]=bloch_siegert(spin_system,ctrl_opers,ctrl_coefs)` | Applies Bloch-Siegert corrections to Cartesian control pulses. Takes control operators and their amplitude arrays and re |
| `kernel/pulses/bruker_write.m` | `bruker_write(X,Y,dt,file_name)` | Saves pulses in Bruker format. The result is a text file with a list of amplitudes and phases, usable in TopSpin. Syntax |
| `kernel/pulses/cartesian2polar.m` | `[r,p,Dr,Dp,Drr,Drp,Dpr,Dpp]=cartesian2polar(x,y,Dx,Dy,Dxx,Dxy,Dyx,Dyy)` | Converts the [RF_x, RF_y] representation of a pulse waveform and the derivatives of any function with respect to those R |
| `kernel/pulses/chirp_pulse.m` | `[Cx,Cy,durs,ints,amps,phis,frqs]=chirp_pulse(npts,dur,bwidth,smp,type)` | Chirp pulse waveform with a sine bell power or a quarter-sine amplitude fade-in and fade-out. Generates unidirectional c |
| `kernel/pulses/grad_pulse.m` | `rho=grad_pulse(spin_system,L,rho,g_amp,s_len,g_dur,s_fac)` | Emulates the effect of a gradient pulse on the sample average density matrix using Edwards formalism. It is assumed that |
| `kernel/pulses/grad_sandw.m` | `rho=grad_sandw(spin_system,L,rho,P,g_amps,s_len,g_durs,s_facs)` | Emulates the effect of a gradient sandwich on the sample average density matrix using Edwards formalism. It is assumed t |
| `kernel/pulses/heterodyne.m` | `[X,Y]=heterodyne(dt,signal,freq)` | Signal heterodyne from wall clock time into the rotating frame using analytic signal demodulation: the negative frequenc |
| `kernel/pulses/isergen.m` | `H=isergen(HL,HM,HR,dt)` | 2nd and 4th order Iserles product quadrature generators for one time propagation step in the case of state-inde- pendent |
| `kernel/pulses/iserstep.m` | `rho_b=iserstep(spin_system,LTM,rho_a,dt)` | Lie-group and Runge-Kutta-Munthe-Kaas solvers for the Lie equa- tion. LG methods are implementations of Equation A.1, wi |
| `kernel/pulses/pmlg5.m` | `phi=pmlg5(n)` | PMLG5 phase sequence as described in the paper by Vinogradova, Madhu and Vega (https://doi.org/10.1016/S0009-2614(99)011 |
| `kernel/pulses/polar2cartesian.m` | `[x,y,Dx,Dy,Dxx,Dxy,Dyx,Dyy]=polar2cartesian(r,p,Dr,Dp,Drr,Drp,Dpr,Dpp)` | Converts [RF_amplitude, RF_phase] representation of a pulse waveform and the derivatives of any function with respect to |
| `kernel/pulses/pulse_demod.m` | `demod_pulse=pulse_demod(time_grid,in_phase,out_phase)` | Interactive demodulation of a complex pulse waveform by a user- specified frequency. Syntax: demod_pulse=pulse_demod(tim |
| `kernel/pulses/pulse_shape.m` | `waveform=pulse_shape(pulse_name,npoints)` | Amplitude envelopes of pulse waveforms. Syntax: waveform=pulse_shape(pulse_name,npoints) Parameters: pulse_name -the nam |
| `kernel/pulses/read_wave.m` | `[A,phi,Cx,Cy,scaling_factor]=read_wave(filename,npoints)` | Reads JCAMP-DX pulse waveform files (a few examples are distri- buted with Spinach, see /kernel/pulses/pk_files). Syntax |
| `kernel/pulses/restrans.m` | `[X,Y,dt]=restrans(X_user,Y_user,dt_user,omega,Q,model,up_factor)` | RLC circuit response calculation -converts a waveform from the ideal shape emitted by the instrument into the shape that |
| `kernel/pulses/rseq_compiler.m` | `[P,T]=rseq_compiler(spin_system,L,Sx,Sy,pulse_phi,pulse_amp,pulse_dur,element_type)` | R sequence compiler. Uses the fact that R-sequences are very repetitive to pre-compile the minimal number of pulse propa |
| `kernel/pulses/rsequence.m` | `[phases,pulse_amp,pulse_dur]=rsequence(n_rotor_periods,n_blocks_per_period,phase_factor,n_cycle_repeats,mas_rate,element_type,supercycle_type)` | R-sequences described in Malcolm Levitt's review: Nomenclature is based on the following notation RN_{n}^{\nu}. Syntax: |
| `kernel/pulses/sawtooth.m` | `waveform=sawtooth(amplitude,frequency,time_grid)` | Returns a saw-tooth waveform. Syntax: waveform=sawtooth(amplitude,frequency,time_grid) Arguments: amplitude - amplitude |
| `kernel/pulses/sech_pulse.m` | `[Cx,Cy,time_grid,amps,phis]=sech_pulse(peak_amp,freq_mod,phase_mod,dur,npts)` | Hyperbolic secant pulse in Cartesian and amplitude-phase representation. Syntax: [Cx,Cy,time_grid,amps,phis]=... sech_pu |
| `kernel/pulses/shaped_pulse_af.m` | `[rho,traj,P]=shaped_pulse_af(spin_system,L0,Lx,Ly,rho,rf_frq_list,rf_amp_list,rf_dur_list,rf_phi,max_rank,method)` | Shaped pulse in amplitude-frequency coordinates using Fokker-Planck formalism (Eqn. 33 in http://dx.doi.org/10.1016/j.jm |
| `kernel/pulses/shaped_pulse_xy.m` | `[rho,traj,P]=shaped_pulse_xy(spin_system,drift,controls,amplitudes,slice_durs,rho,method)` | Shaped pulse function using Cartesian coordinates. Applies a user- specified pulse shape on user-specified operators whi |
| `kernel/pulses/slr_pulse.m` | `[Cx,Cy,durs,amps,phis]=slr_pulse(npts,dur,tbw,flip_angle,pass_rip,stop_rip)` | Shinnar-Le Roux linear-phase selective excitation pulse. Syntax: [Cx,Cy,durs,amps,phis]=slr_pulse(npts,dur,tbw,flip_angl |
| `kernel/pulses/spinal.m` | `phi=spinal(n)` | SPINAL phase sequences as described in the paper by Fung, Khitrin and Ermolaev (https://doi.org/10.1006/jmre.1999.1896). |
| `kernel/pulses/triwave.m` | `waveform=triwave(amplitude,frequency,time_grid)` | Returns a triangular waveform. Syntax: waveform=triwave(amplitude,frequency,time_grid) Arguments: amplitude - amplitude |
| `kernel/pulses/uhrig_times.m` | `time_delays=uhrig_times(T,N)` | Uhrig's UDD decoupling sequence timings. Syntax: time_delays=uhrig_times(T,N) Parameters: T -total duration of the seque |
| `kernel/pulses/vg_pulse.m` | `waveform=vg_pulse(pulse_name,npoints,duration)` | Veshtort-Griffin shaped pulses, generated from tables given in There are good reasons to believe (see Section 2.2 of the |
| `kernel/pulses/wave_basis.m` | `basis_waves=wave_basis(basis_type,n_func,n_points)` | Common basis sets for the expansion of pulse waveforms. Returns the wave- form basis functions as columns of a matrix. S |
| `kernel/reduce.m` | `projectors=reduce(spin_system,L,rho)` | Symmetry and trajectory-level state space reduction. Tries all applicable reduction methods (unless disabled during the |
| `kernel/relaxation.m` | `R=relaxation(spin_system,euler_angles)` | Relaxation superoperator. Syntax: R=relaxation(spin_system,euler_angles) Parameters: euler_angles -three Euler angles (Z |
| `kernel/residual.m` | `spin_system=residual(spin_system)` | Sets up interaction tensors under partial ordering in a liquid crystal with the user-supplied order matrix. All adjustab |
| `kernel/rotframe.m` | `Hr=rotframe(spin_system,H0,H,isotope,order)` | Rotating frame transformation with respect to specified spins to specified order in perturbation theory, using the forma |
| `kernel/spin.m` | `[gamma,multiplicity]=spin(name)` | Database of multiplicities and magnetogyric ratios for sta- ble and long-lived particles, including spin zero. Syntax: [ |
| `kernel/spinlock.m` | `rho=spinlock(spin_system,Lx,Ly,rho,direction)` | Analytical approximation to a spin locking process. This function oblite- rates all spin-spin correlations and all magne |
| `kernel/state.m` | `rho=state(spin_system,states,spins,method)` | Generates Hilbert space density matrices and Liouville space state vectors from their human-readable descriptions. Synta |
| `kernel/states/deut_pair.m` | `[S,T,Q,Tc,Qc]=deut_pair(spin_system,spin_a,spin_b,options)` | All possible states of a spin-1 pair, classified by the total spin into singlet, triplet, and quartet. Syntax: [S,T,Q,Tc |
| `kernel/states/equilibrium.m` | `rho=equilibrium(spin_system,I,Q,euler_angles)` | Returns the thermal equilibrium state at the current temperature. If the anisotropic part and the orientation parameters |
| `kernel/states/four_spin_states.m` | `rho=four_spin_states(spin_system,spins,spin_state)` | Returns user-specified states for a system of four spin-1/2 particles; see also the enclosed Mathematica file. Syntax: r |
| `kernel/states/partner_state.m` | `[A,descr]=partner_state(spin_system,set_spin,partners)` | Partner state expansion; a given state of the specified spins is kroneckered with all combinations of the specified stat |
| `kernel/states/singlet.m` | `S=singlet(spin_system,spin_a,spin_b)` | Returns a two-spin singlet state; both particles must be spin-1/2. Syntax: rho=singlet(spin_system,spin_a,spin_b) Argume |
| `kernel/states/stateinfo.m` | `stateinfo(spin_system,rho,npops)` | Prints the state vector norm and the list of the most populated basis states in the order of decreasing population. Synt |
| `kernel/states/triplet.m` | `[TU,T0,TD]=triplet(spin_system,spin_a,spin_b)` | Returns the components of the two-spin triplet state; both particles must be spin-1/2. Syntax: [Tp,T0,Tm]=triplet(spin_s |
| `kernel/states/unit_state.m` | `rho=unit_state(spin_system)` | Returns the unit state vector or matrix in the current formalism and basis. Syntax: rho=unit_state(spin_system) Paramete |
| `kernel/states/zftrip.m` | `rho=zftrip(spin_system,ZFS,pops,Z,B,idx)` | Projection of the zero-field triplet state with user-specified populations of Cartesian ZFS eigenstates onto the higher- |
| `kernel/steady.m` | `rho=steady(spin_system,P,rho,method)` | Steady state under the repeated action by the same dissi- pative evolution propagator. Syntax: rho=steady(spin_system,P, |
| `kernel/step.m` | `rho=step(spin_system,L,rho,time_step)` | Propagation step function. Computes the action by a matrix exponential without compuing that exponential. Supports one-, |
| `kernel/summaries/summary_basis.m` | `summary_basis(spin_system)` | Prints basis-set state summary for a Spinach system. Syntax: summary_basis(spin_system) Parameters: spin_system -Spinach |
| `kernel/summaries/summary_basis_opts.m` | `summary_basis_opts(spin_system)` | Prints basis-set option summary for a Spinach system. Syntax: summary_basis_opts(spin_system) Parameters: spin_system -S |
| `kernel/summaries/summary_chemistry.m` | `summary_chemistry(spin_system)` | Prints chemical subsystem and exchange summary for a Spinach system. Syntax: summary_chemistry(spin_system) Parameters: |
| `kernel/summaries/summary_coordinates.m` | `summary_coordinates(spin_system,header)` | Prints atomic coordinate summary for a Spinach system. Syntax: summary_coordinates(spin_system,header) Parameters: spin_ |
| `kernel/summaries/summary_couplings.m` | `summary_couplings(spin_system,header)` | Prints spin-spin coupling tensor summary for a Spinach system. Syntax: summary_couplings(spin_system,header) Parameters: |
| `kernel/summaries/summary_mode_coup.m` | `summary_mode_coup(spin_system,header)` | Prints bosonic mode coupling summary for a Spinach system. This covers mode-mode exchange couplings, cross-Kerr coupling |
| `kernel/summaries/summary_mode_mods.m` | `summary_mode_mods(spin_system,header)` | Prints the summary of spin Hamiltonian modulation by bosonic mode coordinates: derivatives of spin-spin coupling tensors |
| `kernel/summaries/summary_modes.m` | `summary_modes(spin_system,header)` | Prints bosonic mode parameter summary for a Spinach system. Syntax: summary_modes(spin_system,header) Parameters: spin_s |
| `kernel/summaries/summary_pbc.m` | `summary_pbc(spin_system,header)` | Prints periodic boundary condition vector summary for a Spinach system. Syntax: summary_pbc(spin_system,header) Paramete |
| `kernel/summaries/summary_rlx_lindblad.m` | `summary_rlx_lindblad(spin_system,header)` | Prints Lindblad relaxation-rate summary for a Spinach system. Syntax: summary_rlx_lindblad(spin_system,header) Parameter |
| `kernel/summaries/summary_rlx_nott.m` | `summary_rlx_nott(spin_system)` | Prints Nottingham DNP relaxation-rate summary for a Spinach system. Syntax: summary_rlx_nott(spin_system) Parameters: sp |
| `kernel/summaries/summary_rlx_t1_t2.m` | `summary_rlx_t1_t2(spin_system,header)` | Prints T1 and T2 relaxation-rate summary for a Spinach system. Syntax: summary_rlx_t1_t2(spin_system,header) Parameters: |
| `kernel/summaries/summary_rlx_weiz.m` | `summary_rlx_weiz(spin_system)` | Prints Weizmann DNP relaxation-rate summary for a Spinach system. Syntax: summary_rlx_weiz(spin_system) Parameters: spin |
| `kernel/summaries/summary_symmetry.m` | `summary_symmetry(spin_system,header)` | Prints permutation-symmetry summary for a Spinach system. Syntax: summary_symmetry(spin_system,header) Parameters: spin_ |
| `kernel/summaries/summary_zeeman.m` | `summary_zeeman(spin_system,header)` | Prints Zeeman interaction tensor summary for a Spinach system. Syntax: summary_zeeman(spin_system,header) Parameters: sp |
| `kernel/thermalize.m` | `R=thermalize(spin_system,R,HLSPS,T,rho_eq,method)` | Modifies the relaxation superoperator to drive the system to the user- specified target state (inhomogeneous master equa |
| `kernel/trajan.m` | `trajan(spin_system,traj,property,time_axis)` | Trajectory analysis function. Plots the time dependence of the densi- ty matrix norm, partitioned into user-specified pr |
| `kernel/trajsimil.m` | `score=trajsimil(spin_system,trajectory_1,trajectory_2,scorefcn)` | Computes trajectory similarity scores. Returns a function representing "similarity" of the two state space trajectories |
| `kernel/utilities/acomm.m` | `C=acomm(A,B)` | A simple shorthand for the anticommutator of two matrices. Syntax: C=acomm(A,B) Parameters: A,B -square matrices Outputs |
| `kernel/utilities/add_spins.m` | `[mult,proj]=add_spins(spin_a,spin_b)` | Reduction of direct products of two su(2) irreps. Syntax: [mult,proj]=add_spins(spin_a,spin_b) Parameters: spin_a -quant |
| `kernel/utilities/adelim.m` | `[L,R]=adelim(spin_system,L,fast_idx,slow_idx)` | Adiabatic elimination in Liouville space, implements Section 6.1 of Kuprov's book. Syntax: [L,R]=adelim(spin_system,L,fa |
| `kernel/utilities/apodisation.m` | `fid=apodisation(spin_system,fid,winfuns,fp_half)` | Performs free induction decay apodisation. Supports free induction decays of any dimension. To satisfy Fourier transform |
| `kernel/utilities/arnoldi.m` | `[V,H]=arnoldi(Op,v0,niter)` | Arnoldi procedure for the creation of an orthonormal Krylov basis from repeated action by an operator on a vector. The p |
| `kernel/utilities/atranspose.m` | `M=atranspose(M)` | Anti-diagonal array transpose. Syntax: M=atranspose(M) Parameters: M -a transposable array Outputs: M -a transposable ar |
| `kernel/utilities/banner.m` | `banner(spin_system,identifier)` | Prints console banners. This is an internal function of the kernel, user calls are discouraged. Syntax: banner(spin_syst |
| `kernel/utilities/binpack.m` | `bins=binpack(box_sizes,bin_size)` | A simple 1D bin packing algorithm. Collects the list of numbers supplied into sublists that sum to the number that is sm |
| `kernel/utilities/blinv.m` | `[Lsq,Dsq]=blinv(A)` | Blicharski's relaxation theory invariants, as given by Equations 20-21 in http://doi.org/10.1515/zna-1972-1012, with an |
| `kernel/utilities/blprod.m` | `[X1_AB,X2_AB]=blprod(A,B)` | Extension of Blicharski's tensor invariants into scalar products of different spin interaction tensors using polarisatio |
| `kernel/utilities/cg_fast.m` | `cg=cg_fast(L,M,L1,M1,L2,M2)` | Clebsch-Gordan coefficient: the coefficient in front of Y(L,M) spheri- cal harmonic in the expansion of the product of Y |
| `kernel/utilities/cheap_norm.m` | `n=cheap_norm(A,t,itmax)` | The cheapest norm for various representations of matrices. CUDA stores matrices by rows, Matlab by columns, and polyadic |
| `kernel/utilities/cheb_coeff.m` | `c=cheb_coeff(f,a,b,n)` | Discrete cosine transform algorithm for Chebyshev expansion coefficients of the user-specified scalar function. Syntax: |
| `kernel/utilities/chemshifts.m` | `[cs_ppm,cs_hz]=chemshifts(spin_system)` | Returns the chemical shifts of every spin in the system relative to the carrier frequency in the current magnet. Syntax: |
| `kernel/utilities/clean_up.m` | `A=clean_up(spin_system,A,nonzero_tol)` | Array clean-up utility. Drops non-zero elements with magnitude below the user-specified tolerance and converts between s |
| `kernel/utilities/clebsch_gordan.m` | `cg=clebsch_gordan(L,M,L1,M1,L2,M2)` | Clebsch-Gordan coefficient: the coefficient in front of Y(L,M) spheri- cal harmonic in the expansion of the product of Y |
| `kernel/utilities/comm.m` | `C=comm(A,B)` | A simple shorthand for the commutator of two matrices. Syntax: C=comm(A,B) Parameters: A,B -square matrices Outputs: C - |
| `kernel/utilities/conmat.m` | `conmatrix=conmat(xyz,r0)` | Molecular connectivity matrix calculator with N*log(N) asymptotic complexity scaling with respect to the num- ber or ato |
| `kernel/utilities/corrfun.m` | `[weights,rates,states]=corrfun(spin_system,n,k,m,p,q)` | Wigner matrix element correlation function under isotropic, axial, and rhombic rotational diffusion. Syntax: [weights,ra |
| `kernel/utilities/cubic_lattice.m` | `[sys,inter]=cubic_lattice(isotope,spacing,n_periods)` | Creates a periodic volume-centered cubic lattice with user- supplied parameters. Syntax: [sys,inter]=cubic_lattice(isoto |
| `kernel/utilities/dfpt.m` | `subgraphs=dfpt(conmatrix,max_sg_size)` | Graph partitioning module. Analyzes the system connectivity graph and creates a list of all connected subgraphs of up to |
| `kernel/utilities/dictum.m` | `spin_system=dictum(spin_system,spins,strength)` | Overrides default assumptions about interaction terms surviving rotating frame transformations. Syntax: spin_system=dict |
| `kernel/utilities/dihedral.m` | `phi=dihedral(A,B,C,D)` | Computes the dihedral angle between vectors specified by the four sets of atomic coordinates. The atoms are assu- med to |
| `kernel/utilities/dilute.m` | `subsystems=dilute(spin_system,isotope,tuples)` | Splits the spin system into several independent subsystems, each containing only one instance of a user specified isotop |
| `kernel/utilities/dipolar.m` | `spin_system=dipolar(spin_system)` | Computes dipolar couplings in the presence or absence of periodic boundary conditions. This is an auxiliary function of |
| `kernel/utilities/dirdiff.m` | `D=dirdiff(spin_system,A,B,T,N)` | Directional derivatives of the matrix exponential. Implements Equation 11 of Najfeld and Havel (https://doi.org/10.1006/ |
| `kernel/utilities/distrib_dim.m` | `A=distrib_dim(A,dim)` | Distributes an array in the user-specified dimension for parallel processing using spmd. Syntax: A=distrib_dim(A,dim) Pa |
| `kernel/utilities/expdrop.m` | `drop=expdrop(from,to,duration,npoints,drop_rate)` | Exponential drop function. Produces an exponential fall-off from a specified value to a specified value with the specifi |
| `kernel/utilities/expmint.m` | `R=expmint(spin_system,A,B,C,T)` | Computes matrix exponential integrals of the following general type: Integrate[expm(-i*A*t)*B*expm(i*C*t),{t,0,T}] Matri |
| `kernel/utilities/expmint2.m` | `I=expmint2(spin_system,A,B,C,D,E,T)` | Computes the nested matrix exponential double integral: Integrate[expm(-i*A*(T-t))*B* Integrate[expm(-i*C*(t-x))*D*expm( |
| `kernel/utilities/fluxonium.m` | `[H,n_op,phi_op]=fluxonium(ec,ej,el,phi_e,nlevels)` | Fluxonium Hamiltonian in a truncated harmonic-oscillator basis. |
| `kernel/utilities/fpl2phan.m` | `phan=fpl2phan(rho,coil,dims)` | Returns the image painted within the Fokker-Planck vector by the user-specified spin state. Syntax: phan=fpl2phan(rho,co |
| `kernel/utilities/fpl2rho.m` | `rho=fpl2rho(rho,dims)` | Returns the average of the spin state vector across the spatial dimensions of the sample. Syntax: rho=fpl2rho(rho,dims) |
| `kernel/utilities/frob_chop.m` | `r=frob_chop(s,tol)` | Truncates SVD decomposition to the user-specified threshold in the Frobenius norm. Syntax: r=frob_chop(s,tol) Parameters |
| `kernel/utilities/g2fplanck.m` | `G=g2fplanck(spin_system,parameters)` | Returns gradient operators within the Fokker-Planck formalism used in the imaging module of Spinach. Syntax: G=g2fplanck |
| `kernel/utilities/get_coupling.m` | `A=get_coupling(spin_system,n,k)` | Extracts the 3x3 coupling tensor between a pair of spins back from the spin_system data structure. Syntax: A=get_couplin |
| `kernel/utilities/gtensorof.m` | `g=gtensorof(spin_system,spin_number)` | Returns the g-tensor of the specified spin at the input orientation. Syntax: g=gtensorof(spin_system,spin_number) Parame |
| `kernel/utilities/hdot.m` | `H=hdot(A,B)` | Hadamard route to Frobenius matrix product. Useful as a replacement for trace(A'*B) because trace(A'*B)=hadm(conj(A),B) |
| `kernel/utilities/herm_spline.m` | `y=herm_spline(f0,df0,f1,df1,x)` | Cubic Hermite spline on [0,1] interval from values and deriva- tives at the interval edges. Syntax: y=herm_spline(f0,df0 |
| `kernel/utilities/human2opspec.m` | `[opspecs,coeffs]=human2opspec(spin_system,operators,spins)` | Converts user-friendly descriptions of spin states and operators into the formal description (opspec) used by Spinach ke |
| `kernel/utilities/hydrodynamics.m` | `[Fx,Fy,Fz]=hydrodynamics(spin_system,parameters)` | A basic hydrodynamics infrastructure provider, returns first derivative operators with respect to the three sample coord |
| `kernel/utilities/idxof.m` | `idx=idxof(sys,label)` | Allows interaction specification by spin label rather than number. Syntax: idx=idxof(sys,label) Parameters: sys -Spinach |
| `kernel/utilities/impound.m` | `answer=impound(varargin)` | This function packages everything it receives into a cell array and returns it back. This is useful for pulling informat |
| `kernel/utilities/intrep.m` | `Hr=intrep(spin_system,H0,H,T,order)` | Interaction representation transformation with respect to a specified Hamiltonian to specified order in perturbation the |
| `kernel/utilities/iselectron.m` | `verdict=iselectron(spin_spec)` | Returns true if the particle is an electron. Syntax: verdict=iselectron(spin_spec) Parameters: spin_spec -a Spinach part |
| `kernel/utilities/iseye.m` | `verdict=iseye(M)` | Returns true for unit matrices. The test is designed to be computationally affordable. Syntax: verdict=iseye(M) Paramete |
| `kernel/utilities/isnucleus.m` | `verdict=isnucleus(spin_spec)` | Returns true if the specification is a nucleus. Syntax: verdict=isnucleus(spin_spec) Parameters: spin_spec -a character |
| `kernel/utilities/isoswap.m` | `[sys,inter]=isoswap(sys,inter,spins,new_iso)` | Makes isotope replacements in the input structures. All interactions are automatically scaled as necessary. Syntax: [sys |
| `kernel/utilities/istraceless.m` | `A=istraceless(M)` | A floating-point precision consistent check for whether a particular matrix is traceless. Syntax: A=istraceless(M) Param |
| `kernel/utilities/isworkernode.m` | `answer=isworkernode()` | Returns true if executed inside a parfor or spmd block. This function is used in the internal decision making of Spinach |
| `kernel/utilities/jacobianest.m` | `[jac,err] = jacobianest(fun,x0)` | Estimate of the Jacobian matrix of a vector valued function of n variables. Syntax: [jac,err] = jacobianest(fun,x0) Para |
| `kernel/utilities/keep_rank.m` | `A=keep_rank(A,nsvk)` | Truncates the singular value decomposition at the specified rank and reassembles the matrix. Syntax: A=keep_rank(A,rank) |
| `kernel/utilities/kill_spin.m` | `spin_system=kill_spin(spin_system,hit_list)` | Removes the specified spins from the spin_system structure and updates it accordingly. Syntax: spin_system=kill_spin(spi |
| `kernel/utilities/killcross.m` | `M=killcross(M,f1idx,f2idx)` | Zeroes the specified rows and columns of a matrix. Syntax: M=killcross(M,f1idx,f2idx) Parameters: M -a matrix f1idx -num |
| `kernel/utilities/killdiag.m` | `spec=killdiag(spec,brush_dim)` | Zeroes out the diagonal of a 2D spectrum using the brush with the specified dimensions. Syntax: spec=killdiag(spec,brush |
| `kernel/utilities/krondelta.m` | `d=krondelta(a,b)` | Kronecker symbol. Syntax: d=krondelta(a,b) Parameters: a -an integer number b -an integer number Outputs: d -a logical n |
| `kernel/utilities/kronm.m` | `x=kronm(Q,x)` | Calculates (Q{1}(x)Q{2}(x)...(x)Q{n})*x without opening Kronecker products. Syntax: x=kronm(Q,x) Parameters: Q - cell ar |
| `kernel/utilities/kronm_new.m` | `M=kronm_new(Q,M)` | Calculates (Q{1}(x)Q{2}(x)...(x)Q{n})*M without opening Kronecker products. Syntax: M=kronm(Q,M) Parameters: Q - cell ar |
| `kernel/utilities/lcurve.m` | `lam_opt=lcurve(lam,err,reg,mode)` | L-curve analysis function. Syntax: lam_opt=lcurve(lam,err,reg,mode) Parameters: lam -row vector of regularisation parame |
| `kernel/utilities/load_vstore.m` | `load_vstore(file_name)` | Loads the current parallel pool ValueStore from a Matlab file. The current store is cleared before the saved keys and va |
| `kernel/utilities/logfactorial.m` | `lf=logfactorial(n)` | Logarithm of the factorial function. Avoids complications with factorials of large numbers overflowing 64-bit numbers. S |
| `kernel/utilities/magpump.m` | `R=magpump(spin_system,R,rho,rate)` | Adds phenomenological pumping terms to the relaxation superoperator to enable approximate simulation of CIDNP, PHIP and |
| `kernel/utilities/md5_hash.m` | `hashstr=md5_hash(A)` | MD5 hash of any Matlab object as a hex string. Identical sparse and full matrices return different hashes. Syntax: hashs |
| `kernel/utilities/merge_inp.m` | `[sys,inter]=merge_inp(sys_parts,inter_parts)` | Merges multiple sys and inter structures into one. Useful for setting up chemical kinetics simulations where the molecul |
| `kernel/utilities/min_int_type.m` | `type=min_int_type(max_val,issigned)` | Minimum integer data type sufficient to store the specified value. Useful in many indexing operati- ons in the Spinach k |
| `kernel/utilities/nearest_spin.m` | `[k,d]=nearest_spin(spin_system,n)` | Returns the index of the nearest spin to the one speci- fied. Only spins for which Cartesian coordinates are available a |
| `kernel/utilities/ngce.m` | `[R,dR]=ngce(spin_system,H0,H1,dt,tau_est,reg)` | Numerical integral route to the Redfield relaxation superopera- tor. Syntax: [R,dR]=ngce(spin_system,H0,H1,dt,tau_est,re |
| `kernel/utilities/nutation_dist.m` | `[freq,distr]=nutation_dist(curve,dt,lambda)` | Nutation frequency distribution from a nutation curve measured with the same coil used for excitation and detection. Syn |
| `kernel/utilities/offsetof.m` | `offs=offsetof(spin_system,idx)` | Returns the isotropic Zeeman offset of the specified spin from the pure magnetogyric ratio frequency in the current magn |
| `kernel/utilities/oscillator.m` | `[H_oscl,X_oscl,xgrid]=oscillator(parameters)` | Harmonic oscillator infrastructure in 1D. Syntax: [H_oscl,X_oscl,xgrid]=oscillator(parameters) Parameters: parameters.fr |
| `kernel/utilities/overwound.m` | `overwound(rho,spc_dim,spn_dim)` | Checks if a Fokker-Planck state vector has any spatial frequencies that its spatial grid is dangerously close to misrepr |
| `kernel/utilities/path_trace.m` | `projectors=path_trace(spin_system,L,rho)` | Liouvillian path tracing. Treats the user-supplied Liouvillian as the adjacency matrix of a graph, computes the weakly c |
| `kernel/utilities/perm_group.m` | `group=perm_group(group_name)` | Permutation group database. Returns complete data for permutation groups. The following group names are available: S2, S |
| `kernel/utilities/phan2fpl.m` | `rho=phan2fpl(phan,rho)` | Projects a spatial intensity distribution into the Fokker-Planck space, using it as the image painted by the the spin st |
| `kernel/utilities/pink_noise.m` | `trajs=pink_noise(amp,dt,npts,ntraj,f_ir,f_uv)` | Synthesises real Gaussian 1/f noise trajectories between specified frequency cutoffs. |
| `kernel/utilities/polinfo.m` | `polinfo(p,level,label)` | Draws an ASCII diagram of a polyadic object. Syntax: polinfo(p) Parameters: p - polyadic object Outputs: an ASCII diagra |
| `kernel/utilities/poolsize.m` | `n=poolsize()` | Returns the current parallel pool size. Syntax: n=poolsize() Parameters: none Outputs: n - number of workers in the curr |
| `kernel/utilities/prune_subgraphs.m` | `subgraphs=prune_subgraphs(subgraphs)` | Removes subgraphs that are contained entirely within other subgraphs. Syntax: subgraphs=prune_subgraphs(subgraphs) Param |
| `kernel/utilities/remncomm.m` | `A=remncomm(A,EvecB,EvalB)` | Removes from the Hermitian operator A the part that does not com- mute with the Hermitian operator B. Syntax: C=remncomm |
| `kernel/utilities/remtrace.m` | `A=remtrace(A)` | Subtracts an appropriate multiple of the unit matrix to make the input matrix traceless. Syntax: A=remtrace(A) Parameter |
| `kernel/utilities/repcols.m` | `B=repcols(A,col_nums,rep_counts)` | Replicates specified columns of a matrix or cell array a specified number of times. Syntax: B=repcols(A,col_nums,rep_cou |
| `kernel/utilities/report.m` | `report(spin_system,report_string)` | Writes a log message to the console or an ACSII file. The message includes the call stack of the function that produced |
| `kernel/utilities/reprows.m` | `B=reprows(A,row_nums,rep_counts)` | Replicates specified rows of a matrix or cell array a specified number of times. Syntax: B=reprows(A,row_nums,rep_counts |
| `kernel/utilities/rlx_modes.m` | `R=rlx_modes(spin_system)` | Bosonic mode dissipation superoperator. Builds thermalised GKSL dissipators for the amplitude damping and the pure depha |
| `kernel/utilities/rlx_phonon.m` | `R=rlx_phonon(spin_system,H,X,I0,alpha,T,form)` | Spin-phonon relaxation in the generalised Lindblad form of Saito, Miyashita, and De Raedt (Phys. Rev. B 60, 14553 (1999)), |
| `kernel/utilities/rlx_scalar.m` | `R=rlx_scalar(spin_system,H0,H1,tau_c_array)` | Scalar relaxation superoperator using Redfield theory. Syntax: R=rlx_scalar(spin_system,H0,H1,tau_c_array) Parameters: H |
| `kernel/utilities/rlx_split.m` | `[R1,R2,Rm]=rlx_split(spin_system,R)` | Splits a relaxation superoperator into longitudinal, trans- verse and mixed components. Syntax: [R1,R2,Rm]=rlx_split(spi |
| `kernel/utilities/rlx_t1_t2.m` | `[R1Op,R2Op]=rlx_t1_t2(spin_system,euler_angles)` | Extended T1/T2 relaxation model returning the relaxation super- operators separately for the longitudinal and the transv |
| `kernel/utilities/rocomm.m` | `C=rocomm(A)` | Right-ordered nested commutator [[[[A{1},A{2}],A{3}],A{4}],...] built from the user-supplied matrices. Syntax: C=rocomm( |
| `kernel/utilities/rotor_stack.m` | `[L,rotor_phases]=rotor_stack(spin_system,parameters,assumptions)` | Returns a rotor stack of Liouvillians or Hamiltonians. The stack is needed for the traditional style calculation of MAS |
| `kernel/utilities/rspert.m` | `[Ep,Vp]=rspert(E0,H1,order)` | Rayleigh-Schrodinger perturbation theory to arbitrary order, Eqs 2.21-2.23 from Stefan Stoll's PhD thesis, with the typo |
| `kernel/utilities/rspt_eig.m` | `[E,V,dE,T,LP]=rspt_eig(spin_system,parameters,Hz,Hc,Hmw,B)` | Eigensystem of sparse Hamiltonians to user-specified order in RSPT with careful handling of diagonal dominance and an op |
| `kernel/utilities/rwalk.m` | `eulers=rwalk(npts,tau_c,dt)` | Random walk on SO(3), isotropic rotational diffusion. Syntax: eulers=rwalk(npts,tau_c,dt) Parameters: npts -number of po |
| `kernel/utilities/save_vstore.m` | `save_vstore(file_name)` | Saves the current parallel pool ValueStore into a Matlab file. The snapshot contains keys and values only; callback func |
| `kernel/utilities/scomponents.m` | `sci=scomponents(A)` | Strongly connected components of a graph, David Gleich's imple- mentation of Tarjan's algorithm: Syntax: sci=scomponents |
| `kernel/utilities/sec2kite.m` | `R=sec2kite(spin_system,R)` | Converts a secular relaxation superoperator into the Redfield kite form by dropping all non-longitudinal cross-relaxatio |
| `kernel/utilities/shift_iso.m` | `tensors=shift_iso(tensors,spin_numbers,new_iso)` | Replaces the isotropic parts of interaction tensors with user- supplied values. This is useful for correcting DFT calcul |
| `kernel/utilities/sim2liouv.m` | `[spin_system,parameters,H,R,K]=sim2liouv(spin_system,parameters,H,R,K)` | Moves a zeeman-hilb simulation context into Liouville space. When the formalism specified in the spin system object is ' |
| `kernel/utilities/sinkhole.m` | `L=sinkhole(spin_system,L,states)` | Turns the specified states into sinkholes --any population reaching them will be summed up and stored forever in a froze |
| `kernel/utilities/snormpdf.m` | `p=snormpdf(x,mu,sigma,alpha)` | Azzalini's skew normal distribution. Syntax: p=snormpdf(x,mu,sigma,alpha) Parameters: x -an array of real numbers mu -ex |
| `kernel/utilities/sorensen.m` | `b=sorensen(rho_init,rho_targ)` | Sorensen bound for the maximum transfer efficiency between two states under arbitrary control operators. Equation 186 fr |
| `kernel/utilities/sp_block_diag.m` | `S=sp_block_diag(varargin)` | Sparse block diagonal matrix from a stack of matrix blocks. Syntax: S=sp_block_diag(A) S=sp_block_diag(A,B,C,...) Parame |
| `kernel/utilities/sparse2csr.m` | `[row_ptr,col_idx]=sparse2csr(A)` | Computes a partial compressed row storage (CSR) transformation for a given Matlab sparse matrix. Adapted from the code w |
| `kernel/utilities/spden.m` | `J=spden(L,D,omega)` | Lorentzian spectral density function for rotational diffusion at the user-specified frequency. Syntax: J=spden(L,D,omega |
| `kernel/utilities/spher_harmon.m` | `Y=spher_harmon(l,m,theta,phi)` | Spherical harmonics. Syntax: Y=spher_harmon(l,m,theta,phi) Parameters: l -L quantum number m -M quantum number theta -an |
| `kernel/utilities/sphten2zeeman.m` | `P=sphten2zeeman(spin_system)` | Returns a matrix that converts state vectors written in the spherical tensor basis set used by Spinach into state vector |
| `kernel/utilities/stitch.m` | `fid=stitch(spin_system,L,rho_stack,coil_stack,mec_oper,mec_time,t1,t2,t3,tdir)` | Stitching function for bidirectionally propagated 3D NMR pulse sequences. Propagate your initial condition forward to so |
| `kernel/utilities/svd_shrink.m` | `[vec,cov]=svd_shrink(spin_system,rho,tol)` | Generates sets of vector-covector pairs for the parallel implementation of the time propagation algorithm described in [ |
| `kernel/utilities/swizzle.m` | `tuples=swizzle(index_arrays)` | Flattens out nested index lists and outputs them as an array of tuples in random order. This is useful for flattening ne |
| `kernel/utilities/symmetry.m` | `spin_system=symmetry(spin_system,bas)` | Permutation symmetry treatment. Compiles character tables of composite symmetry groups, builds the permutation table for |
| `kernel/utilities/tikhoind.m` | `[x,err,reg]=tikhoind(K,D,y,lam)` | Analytical Tikhonov regularised solution to K*x=y without any constraints (sign-indefinite output). Syntax: [x,err,reg]= |
| `kernel/utilities/tikhol1n.m` | `[x,err,reg]=tikhol1n(A,y,nnzt)` | L1 norm Tikhonov regularised solver for A*x=y where A is an ill-conditioned matrix. The error functional is norm(A*x-y,2 |
| `kernel/utilities/tikhonov.m` | `[x,err,reg]=tikhonov(K,D,KtK,DtD,H,y,lambda)` | Tikhonov regularised solution to K*x=y with a positivity const- raint on x using regularised Newton-Raphson method. Synt |
| `kernel/utilities/tolerances.m` | `[spin_system,sys]=tolerances(spin_system,sys)` | Tolerances and fundamental constants. Sets various accuracy cut-offs, constants and tolerances used by Spinach kernel. S |
| `kernel/utilities/transfermat.m` | `T=transfermat(amp_inps,amp_outs)` | Transfer matrix calculation for linear filters. Syntax: T=transfermat(amp_inps,amp_outs) Parameters: amp_inps -a matrix |
| `kernel/utilities/unihash.m` | `A=unihash(A)` | Hash table based stable duplicate row eliminator, for use with large sparse matrices where Matlab's unique(...,'rows') i |
| `kernel/utilities/v2fplanck.m` | `F=v2fplanck(spin_system,parameters)` | Translates a stationary 3D velocity field and a diffusion tensor field into a Fokker-Planck evolution generator. Syntax: |
| `kernel/utilities/validate_sym.m` | `validate_sym(spin_system,bas)` | Extended validation of user-declared permutation symmetry. Confirms that the Zeeman, coupling, and giant-spin interactio |
| `kernel/utilities/vvpert.m` | `[Ep,G]=vvpert(E0,H1,order)` | Van Vleck perturbation theory, following Shavitt and Redmon, but excluding the quasi-degenerate split. Syntax: [Ep,G]=vv |
| `kernel/utilities/which_subst.m` | `subst=which_subst(spin_system,spins)` | Finds out which substance hosts the specified spins; throws an error if there is more than one. Syntax: subst=which_subs |
| `kernel/utilities/wigner.m` | `D=wigner(l,alp,bet,gam)` | Wigner D matrices, defined as (Brink & Satchler, Eq 2.13): D=expm(-1i*Lz*alp)*expm(-1i*Ly*bet)*expm(-1i*Lz*gam); where L |
| `kernel/utilities/wigner_3j.m` | `w=wigner_3j(j1,m1,j2,m2,j3,m3)` | Calculates Wigner 3j-symbols. Syntax: w=wigner_3j(j1,m1,j2,m2,j3,m3) If physically inadmissible indices are supplied, a |
| `kernel/utilities/wigner_6j.m` | `w=wigner_6j(j1,j2,j3,j4,j5,j6)` | Wigner 6j-symbols. Syntax: w=wigner_6j(j1,j2,j3,j4,j5,j6) If physically inadmissible indices are supplied, a zero is ret |
| `kernel/utilities/wigner_fock.m` | `W=wigner_fock(rho,alpha)` | Wigner function of a bosonic mode state given as a density matrix in a truncated Fock basis, evaluated at the specified |
| `kernel/utilities/xyz2dd.m` | `[d,alp,bet,gam,M]=xyz2dd(r1,r2,isotope1,isotope2)` | Converts coordinate specification of the dipolar interaction into the dipolar interaction constant, three Euler angles, |
| `kernel/utilities/xyz2hfc.m` | `A=xyz2hfc(exyz,nxyz,isotope)` | Converts point electron and nuclear coordinates into a hyper- fine interaction tensor. Syntax: A=xyz2hfc(exyz,nxyz,isoto |
| `kernel/utilities/xyz2pd.m` | `density=xyz2pd(coords,x_range,y_range,z_range,x_npts, y_npts, z_npts)` | Probability density estimation for a three-dimensional Cartesian point cloud on a user-specified regular grid. Syntax: d |
| `kernel/utilities/zte.m` | `projector=zte(spin_system,L,rho,nstates)` | Zero track elimination function. Inspects the first few steps in the system trajectory and drops the states that did not |
