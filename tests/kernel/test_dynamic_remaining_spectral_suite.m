% Tests remaining spectral, symmetry, and Fokker-Planck utilities. Syntax:
%
%                    result=test_dynamic_remaining_spectral_suite()
%
% Outputs:
%
%     result  - regression test result with explanatory messages
%
% The test checks field-swept eigensystem helpers, rotor-stack assembly,
% permutation-symmetry projectors, and one-dimensional Fokker-Planck
% operators against compact analytical references.
%
% ilya.kuprov@weizmann.ac.il

function result=test_dynamic_remaining_spectral_suite()

% Announce the test target
fprintf('TESTING: Remaining spectral and spatial utilities\n');

% State the utility target of the test
result=new_test_result('kernel/dynamic_remaining_spectral_suite',...
                       'Remaining spectral and spatial utilities',...
                       'Small spectral, symmetry, and Fokker-Planck helpers must match analytical matrix references.');

% Check exact diagonalisation path in rspt_eig against eig()
spin_system=local_minimal_system('zeeman-hilb',2);
parameters.rspt_order=Inf;
parameters.rho0=diag([0.9 0.1]);
Hz=diag([1 3]);
Hc=[0 0.1;0.1 0];
Hmw=[0 1;1 0];
field=2;
[E_obs,~,dE_obs,T_obs,LP_obs]=rspt_eig(spin_system,parameters,Hz,Hc,Hmw,field);
[V_ref,E_ref]=eig(field*Hz+Hc,'vector');
[E_ref,sort_idx]=sort(E_ref,'ascend');
V_ref=V_ref(:,sort_idx);
dE_ref=real(diag(V_ref'*Hz*V_ref));
T_ref=abs(V_ref'*Hmw*V_ref).^2;
LP_ref=real(diag(V_ref'*parameters.rho0*V_ref));
result=test_close(result,'rspt_eig exact energies',E_obs,E_ref,1e-14,1e-14,...
                  'the infinite-order branch must agree with dense Hermitian diagonalisation');
result=test_close(result,'rspt_eig Hellmann-Feynman derivatives',dE_obs,dE_ref,1e-14,1e-14,...
                  'dE/dB must be the Zeeman operator expectation value in each eigenstate');
result=test_close(result,'rspt_eig transition moments',T_obs,T_ref,1e-14,1e-14,...
                  'transition moments are squared microwave-operator matrix elements');
result=test_close(result,'rspt_eig level populations',LP_obs,LP_ref,1e-14,1e-14,...
                  'level populations are equilibrium density-matrix expectations');

% Check Liouville-space resonance-field extraction for a diagonal pencil
spin_system=local_minimal_system('sphten-liouv',2);
parameters=struct();
parameters.mw_freq=100;
parameters.window=[9 11];
parameters.pp_tol=1e-9;
parameters.tm_tol=0.25;
parameters.fwhm=0.01;
parameters.orientation=[0 0 0];
Iz=2*pi*diag([10 20]);
Ic=zeros(2);
Qz=local_tiny_rank_one(2);
Qc=local_tiny_rank_one(2);
Hmw=[1;0];
[tf,tm,tw,pd]=eigenfields(spin_system,parameters,Iz,Qz,Ic,Qc,Hmw);
result=test_close(result,'eigenfields transition field',tf,10,1e-8,1e-12,...
                  'a 100 Hz transition under a 10 Hz/T Liouville pencil occurs at 10 T');
result=test_close(result,'eigenfields transition moment',tm,1,1e-14,1e-14,...
                  'the selected normalised dyadic has unit microwave transition moment');
result=test_close(result,'eigenfields transition width',tw,parameters.fwhm,1e-14,1e-14,...
                  'Liouville-space eigenfields return the requested phenomenological FWHM');
result=test_close(result,'eigenfields population difference',pd,1,1e-14,1e-14,...
                  'Liouville-space population differences are currently unit placeholders');

% Check one-dimensional gradient operator construction in Fokker-Planck space
spin_system=local_created_system('zeeman-hilb',1);
fp_par.dims=2;
fp_par.npts=3;
G=g2fplanck(spin_system,fp_par);
H=hamiltonian(assume(spin_system,'labframe','zeeman'))/spin_system.inter.magnet;
G_ref=kron(spdiags([-1;0;1],0,3,3),H);
result=test_close(result,'g2fplanck one-dimensional gradient',inflate(G{1}),G_ref,1e-12,1e-12,...
                  'the 1D gradient operator is the centred spatial coordinate Kronecker H per tesla');
result=test_true(result,'g2fplanck empty inactive dimensions',isempty(G{2})&&isempty(G{3}),...
                 'inactive spatial dimensions must return empty gradient operators');

% Check one-dimensional velocity generator against the finite-difference derivative
spin_system=local_minimal_system('sphten-liouv',1);
sp_par.dims=1;
sp_par.npts=10;
sp_par.deriv={'period',3};
sp_par.u=0.2*ones(10,1);
F_obs=v2fplanck(spin_system,sp_par);
Dx=fdmat(sp_par.npts,sp_par.deriv{2},1)/(sp_par.dims/sp_par.npts);
Fx=-1i*Dx;
F_ref=spdiags(Fx*sp_par.u(:),0,10,10)+spdiags(sp_par.u(:),0,10,10)*Fx;
F_ref=clean_up(spin_system,F_ref,spin_system.tols.liouv_zero);
result=test_close(result,'v2fplanck one-dimensional flow',F_obs,F_ref,1e-12,1e-12,...
                  'a uniform 1D velocity field must assemble the expected finite-difference flow generator');

% Check rotor-stack rank-zero assembly against direct Hamiltonian generation
spin_system=local_created_system('sphten-liouv',14.1);
mas_par.axis=[0 0 1];
mas_par.offset=0;
mas_par.spins={'1H'};
mas_par.max_rank=0;
mas_par.rframes={};
mas_par.orientation=[0 0 0];
mas_par.masframe='rotor';
local_ensure_pool();
[L_stack,rotor_phases]=rotor_stack(spin_system,mas_par,'nmr');
[H_direct,~]=hamiltonian(assume(spin_system,'nmr'));
H_direct=frqoffset(spin_system,H_direct,mas_par);
H_direct=clean_up(spin_system,H_direct,spin_system.tols.liouv_zero);
result=test_close(result,'rotor_stack rank-zero Hamiltonian',L_stack{1},H_direct,1e-12,1e-12,...
                  'with max_rank zero and no anisotropy, rotor_stack must return the direct Hamiltonian');
result=test_close(result,'rotor_stack phase grid',rotor_phases,0,1e-14,1e-14,...
                  'a rank-zero rotor stack has a single zero rotor phase');

% Check S2 fully symmetric projector construction on a four-state product basis
spin_system=local_minimal_system('sphten-liouv',4);
spin_system.comp.nspins=2;
spin_system.bas.basis=[0 0;1 0;0 1;1 1];
bas.sym_group={'S2'};
bas.sym_spins={[1 2]};
bas.sym_a1g_only=true;
spin_system=symmetry(spin_system,bas);
projector=spin_system.bas.irrep.projector;
result=test_close(result,'symmetry projector orthonormality',projector'*projector,eye(3),...
                  1e-14,1e-14,...
                  'A1g orbit projectors should be orthonormal after orbit normalisation');
result=test_close(result,'symmetry mixed orbit',abs(projector(:,2)),[0;1;1;0]/sqrt(2),...
                  1e-14,1e-14,...
                  'the mixed S2 orbit must symmetrise the two exchanged basis states');
result=test_close(result,'symmetry irrep dimension',spin_system.bas.irrep.dimension,3,...
                  1e-14,1e-14,...
                  'the four-state basis has three S2 orbit representatives in the A1g irrep');

end


function spin_system=local_created_system(formalism,magnet)

% Build a quiet one-proton Spinach object in the requested formalism
sys.magnet=magnet;
sys.isotopes={'1H'};
inter.zeeman.scalar={0};
bas.formalism=formalism;
bas.approximation='none';
if strcmp(formalism,'sphten-liouv')
    bas.projections=+1;
end
spin_system=test_spin_system(sys,inter,bas);

end


function spin_system=local_minimal_system(formalism,dim)

% Create a quiet minimal descriptor for matrix-only helper calls
spin_system.sys.output='hush';
spin_system.sys.enable={};
spin_system.sys.disable={};
spin_system.bas.formalism=formalism;
spin_system.bas.basis=zeros(dim,1);
spin_system.comp.mults=2;
spin_system.tols.liouv_zero=1e-12;
spin_system.tols.prop_chop=1e-12;
spin_system.tols.dense_matrix=0.5;
spin_system.tols.small_matrix=10;

end


function Q=local_tiny_rank_one(dim)

% Create a rank-one rotational-basis cell that evaluates to a tiny zero surrogate
Q=cell(1,1);
Q{1}=cell(3,3);
Q{1}{2,2}=sparse(1,1,realmin,dim,dim);

end


function local_ensure_pool()

% Start a one-worker process pool for compact parfor utilities
current_pool=gcp('nocreate');
if isempty(current_pool)
    parpool('Processes',1);
end

end

