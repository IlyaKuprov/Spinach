% Tests state-constructor helper functions. Syntax:
%
%                    result=test_state_constructor_suite()
%
% Outputs:
%
%     result  - regression test result with explanatory messages
%
% The test checks unit, thermal, singlet/triplet, partner, deuteron-pair,
% four-spin, and zero-field triplet state constructors using projector and
% normalisation identities.
%
% ilya.kuprov@weizmann.ac.il

function result=test_state_constructor_suite()

% Announce the test target
fprintf('TESTING: State-constructor functions\n');

% State the state-constructor target of the test
result=new_test_result('kernel/state_constructor_suite',...
                       'State-constructor functions',...
                       'state constructors must produce correctly normalised physical density objects.');

% One-spin Hilbert and Liouville unit states have known forms
sys.magnet=0;
sys.isotopes={'1H'};
inter.zeeman.scalar={0};
inter.temperature=300;
bas.formalism='zeeman-hilb'; bas.approximation='none';
spin_h=test_spin_system(sys,inter,bas);
result=test_close(result,'unit_state zeeman-hilb',unit_state(spin_h),speye(2),1e-15,1e-15,...
                  'Hilbert-space unit state is the sparse unit density matrix');
result=test_close(result,'equilibrium zero Hamiltonian',equilibrium(spin_h,zeros(2)),speye(2)/2,1e-14,1e-14,...
                  'a zero Hamiltonian has the maximally mixed thermal state at finite temperature');
bas.formalism='zeeman-liouv';
spin_l=test_spin_system(sys,inter,bas);
unit_zeeman=speye(2); unit_zeeman=unit_zeeman(:)/sqrt(2);
result=test_close(result,'unit_state zeeman-liouv',unit_state(spin_l),unit_zeeman,1e-15,1e-15,...
                  'Zeeman Liouville unit state is the Frobenius-normalised stretched unit matrix');
bas.formalism='sphten-liouv';
spin_s=test_spin_system(sys,inter,bas);
unit_s=unit_state(spin_s);
result=test_close(result,'unit_state sphten-liouv first component',unit_s,sparse(1,1,1,size(spin_s.bas.basis,1),1),1e-15,1e-15,...
                  'spherical-tensor Liouville unit state is the T(0,0) population basis vector');
stateinfo(spin_s,unit_s,1);
result=test_true(result,'stateinfo sphten-liouv smoke',true,...
                 'stateinfo must accept a valid spherical-tensor state vector and population count');

% Two-spin singlet and triplet constructors must form four orthogonal projectors summing to identity
sys2.magnet=0;
sys2.isotopes={'1H','1H'};
inter2.zeeman.scalar={0,0};
bas2.formalism='zeeman-hilb'; bas2.approximation='none';
spin2=test_spin_system(sys2,inter2,bas2);
S=singlet(spin2,1,2);
[TU,T0,TD]=triplet(spin2,1,2);
result=test_close(result,'singlet trace',trace(S),1,1e-14,1e-14,...
                  'a pure two-spin singlet density matrix has unit trace');
result=test_close(result,'triplet traces',[trace(TU) trace(T0) trace(TD)],[1 1 1],1e-14,1e-14,...
                  'each triplet projection is a unit-trace density matrix');
result=test_close(result,'singlet triplet completeness',S+TU+T0+TD,speye(4),1e-14,1e-14,...
                  'singlet plus three triplet projectors resolve the two-spin identity');
result=test_close(result,'singlet projector idempotence',S*S,S,1e-14,1e-14,...
                  'the singlet density matrix is a rank-one projector');
result=test_close(result,'triplet zero orthogonality',trace(S*T0),0,1e-14,1e-14,...
                  'singlet and triplet subspaces are orthogonal');

% partner_state must enumerate all requested partner-state combinations
sys3.magnet=0;
sys3.isotopes={'1H','1H','1H'};
inter3.zeeman.scalar={0,0,0};
spin3=test_spin_system(sys3,inter3,bas2);
[A,descr]=partner_state(spin3,{{'L+',2}},{{{'E','Lz'},[1 3]}});
result=test_close(result,'partner_state combination count',numel(A),4,0,0,...
                  'two binary partner spins generate four product-state combinations');
for n=1:numel(A)
    result=test_close(result,['partner_state descriptor ' int2str(n)],A{n},state(spin3,descr{n},{1,2,3}),1e-14,1e-14,...
                      'each partner_state descriptor must reproduce the corresponding Spinach state');
end

% Four-spin singlet-singlet state must match the explicit product of two two-spin singlets
sys4.magnet=0;
sys4.isotopes={'1H','1H','1H','1H'};
inter4.zeeman.scalar={0,0,0,0};
spin4=test_spin_system(sys4,inter4,bas2);
spin2a=test_spin_system(sys2,inter2,bas2);
SS_ref=kron(singlet(spin2a,1,2),singlet(spin2a,1,2));
SS=four_spin_states(spin4,1:4,'S(x)S');
result=test_close(result,'four_spin_states SxS',SS,SS_ref,1e-14,1e-14,...
                  'the four-spin S(x)S state is the direct product of two spin-half singlet projectors');

% Spin-1 pair states must resolve the nine-dimensional Hilbert identity
sysd.magnet=0;
sysd.isotopes={'2H','2H'};
interd.zeeman.scalar={0,0};
spind=test_spin_system(sysd,interd,bas2);
[Sd,Td,Qd]=deut_pair(spind,1,2);
projector_sum=Sd;
for n=1:numel(Td), projector_sum=projector_sum+Td{n}; end
for n=1:numel(Qd), projector_sum=projector_sum+Qd{n}; end
result=test_close(result,'deut_pair completeness',projector_sum,speye(9),1e-12,1e-12,...
                  'spin-1 singlet, triplet, and quintet projectors resolve the pair Hilbert space');
result=test_close(result,'deut_pair traces',[trace(Sd) cellfun(@trace,Td) cellfun(@trace,Qd)],[1 ones(1,3) ones(1,5)],1e-12,1e-12,...
                  'all deuteron-pair population states are unit-trace projectors');

% Zero-field triplet projection must preserve density-matrix trace and Hermiticity
syse.magnet=0;
syse.isotopes={'E3'};
intere.zeeman.scalar={0};
spine=test_spin_system(syse,intere,bas2);
ZFS=diag([-1 0 1])*1e6;
Z=eye(3)*28e9;
rho_zf=zftrip(spine,ZFS,[0.2 0.3 0.5],Z,0.01,1);
result=test_close(result,'zftrip trace',trace(rho_zf),1,1e-12,1e-12,...
                  'zero-field triplet population projection must preserve unit trace');
result=test_close(result,'zftrip Hermiticity',rho_zf,rho_zf',1e-12,1e-12,...
                  'projected zero-field triplet density matrix must be Hermitian');

end
