% Tests deterministic physical formula utility helpers. Syntax:
%
%                    result=test_dynamic_physics_formulae_suite()
%
% Outputs:
%
%     result  - regression test result with explanatory messages
%
% The test checks spin addition, point-dipole tensors, hyperfine tensors,
% exponential drops, skew-normal densities, oscillator grids, hydrodynamic
% derivative construction, and spherical-tensor projection metadata.
%
% ilya.kuprov@weizmann.ac.il

function result=test_dynamic_physics_formulae_suite()

% State the physical formula target of the test
result=new_test_result('kernel/dynamic_physics_formulae_suite',...
                       'Physical formula utilities',...
                       'closed-form physical helper formulae must match their analytic definitions on small systems.');

% Check Clebsch-Gordan reduction of two spin-half irreps
[mults,proj]=add_spins(1/2,1/2);
P=[proj{1} proj{2}];
result=test_true(result,'add_spins multiplicities',isequal(mults,[1 3]),...
                 'two spin-half irreps reduce into one singlet and one triplet irrep');
result=test_close(result,'add_spins projector unity',P'*P,eye(4),1e-12,1e-12,...
                  'canonical projectors returned by add_spins must be orthonormal');
result=test_close(result,'add_spins complete basis',P*P',eye(4),1e-12,1e-12,...
                  'singlet and triplet projectors must span the four-dimensional product space');

% Check point-dipole coupling for a one-Angstrom z-axis displacement
hbar=1.054571628e-34;
mu0=4*pi*1e-7;
[d,alp,bet,gam,M]=xyz2dd([0 0 0],[0 0 1],'1H','13C');
d_ref=spin('1H')*spin('13C')*hbar*mu0/(4*pi*(1e-10)^3);
result=test_close(result,'xyz2dd coupling constant',d,d_ref,1e-6,1e-12,...
                  'the point-dipole constant is gamma1*gamma2*hbar*mu0/(4*pi*r^3)');
result=test_close(result,'xyz2dd euler angles',[alp bet gam],[0 0 0],1e-14,1e-14,...
                  'a z-axis dipolar tensor is already in its principal-axis frame');
result=test_close(result,'xyz2dd tensor',M,d_ref*diag([1 1 -2]),1e-6,1e-12,...
                  'a z-axis point dipole has traceless principal values d, d, and -2d');

% Check point electron-nucleus hyperfine tensor for a z-axis displacement
A=xyz2hfc([0 0 0],[0 0 1],'1H',1);
C=1e4*spin('1H')*hbar*mu0/(4*pi*(1e-10)^3);
result=test_close(result,'xyz2hfc point tensor',A,C*diag([-1 -1 2]),1e-4,1e-12,...
                  'a z-axis point electron gives a Gauss hyperfine tensor proportional to diag(-1,-1,2)');

% Check exponential drop values at exact quartering points
rate=log(4);
drop=expdrop(5,1,2,3,rate);
result=test_close(result,'expdrop exact points',drop,[5 9/5 1],1e-14,1e-14,...
                  'with rate log(4) over two seconds the middle exponential factor is one quarter');

% Check skew-normal density reduces to the normal density at zero skew
x=[-1 0 1];
p=snormpdf(x,0,1,0);
p_ref=exp(-(x.^2)/2)/sqrt(2*pi);
result=test_close(result,'snormpdf zero skew',p,p_ref,1e-15,1e-15,...
                  'Azzalini skew-normal density is the ordinary normal density when alpha is zero');

% Check oscillator coordinate grid and coordinate operator
parameters.frc_cnst=4;
parameters.par_mass=2;
parameters.grv_cnst=0;
parameters.n_points=11;
parameters.box_size=10;
[H_oscl,X_oscl,xgrid]=oscillator(parameters);
result=test_close(result,'oscillator grid',xgrid,(-5:5)',1e-15,1e-15,...
                  'oscillator grid points must span the box uniformly and symmetrically');
result=test_close(result,'oscillator coordinate operator',diag(X_oscl),xgrid,1e-15,1e-15,...
                  'the oscillator coordinate operator must have the grid on its diagonal');
result=test_close(result,'oscillator hermiticity',H_oscl-H_oscl',zeros(11),1e-15,1e-15,...
                  'the finite-difference oscillator Hamiltonian must be Hermitian');

% Check one-dimensional hydrodynamic derivative construction
hydro_spin_system.sys.enable={'polyadic'};
hydro_params.dims=10;
hydro_params.npts=10;
hydro_params.deriv={'period',3};
[Fx,Fy,Fz]=hydrodynamics(hydro_spin_system,hydro_params);
Dx=fdmat(10,3,1)/(hydro_params.dims/hydro_params.npts);
result=test_close(result,'hydrodynamics one-dimensional Fx',inflate(Fx),-1i*Dx,1e-15,1e-15,...
                  'one-dimensional hydrodynamics must wrap the finite-difference derivative in a polyadic');
result=test_true(result,'hydrodynamics empty transverse axes',isempty(Fy)&&isempty(Fz),...
                 'unused y and z hydrodynamic derivative operators must be empty in one dimension');

% Check spherical-tensor to Zeeman projection metadata on one spin-half
spin_system.comp.mults=2;
spin_system.bas.formalism='sphten-liouv';
spin_system.bas.basis=(0:3)';
P=sphten2zeeman(spin_system);
result=test_true(result,'sphten2zeeman dimensions',isequal(size(P),[4 4])&&rank(full(P))==4,...
                 'one spin-half has four spherical-tensor basis states and four Zeeman-Liouville states');

end


