% Tests dynamic state projection helper paths. Syntax:
%
%              result=test_dynamic_state_projection_suite()
%
% Outputs:
%
%     result  - regression test result with explanatory messages
%
% The test checks deuteron-pair coherences, dephased population stationarity,
% captured stateinfo() output, and isotropic zero-field triplet projection.
%
% ilya.kuprov@weizmann.ac.il

function result=test_dynamic_state_projection_suite()

% State the state-projection target of the test
result=new_test_result('kernel/dynamic_state_projection_suite',...
                       'Dynamic state projection helpers',...
                       'state projection helpers must preserve Hermitian adjoint pairs, stationarity, and normalisation.');

% Build a two-deuteron Hilbert-space system
sys.magnet=0;
sys.isotopes={'2H','2H'};
inter.zeeman.scalar={0,0};
inter.temperature=300;
bas.formalism='zeeman-hilb';
bas.approximation='none';
spin_d=test_spin_system(sys,inter,bas);

% Request all deuteron-pair coherences
[~,~,~,Tc,Qc]=deut_pair(spin_d,1,2);

% Check adjoint pairing of triplet coherences
result=test_close(result,'triplet coherence adjoint 1',Tc{1},Tc{3}',1e-12,1e-12,...
                  'T0->T- and T-->T0 coherences must be Hermitian adjoints');
result=test_close(result,'triplet coherence adjoint 2',Tc{2},Tc{4}',1e-12,1e-12,...
                  'T+->T0 and T0->T+ coherences must be Hermitian adjoints');

% Check adjoint pairing of quintet coherences
result=test_close(result,'quintet coherence adjoint 1',Qc{1},Qc{5}',1e-12,1e-12,...
                  'Q-->Q-- and Q--->Q- coherences must be Hermitian adjoints');
result=test_close(result,'quintet coherence adjoint 2',Qc{2},Qc{6}',1e-12,1e-12,...
                  'Q0->Q- and Q-->Q0 coherences must be Hermitian adjoints');
result=test_close(result,'quintet coherence adjoint 3',Qc{3},Qc{7}',1e-12,1e-12,...
                  'Q+->Q0 and Q0->Q+ coherences must be Hermitian adjoints');
result=test_close(result,'quintet coherence adjoint 4',Qc{4},Qc{8}',1e-12,1e-12,...
                  'Q++->Q+ and Q+->Q++ coherences must be Hermitian adjoints');

% Request dephased deuteron-pair populations
options.dephasing=1;
[Sd,Td,Qd]=deut_pair(spin_d,1,2,options);

% Build the longitudinal two-spin Hamiltonian used by the dephasing claim
Hz=operator(spin_d,'Lz',1)+operator(spin_d,'Lz',2);

% Check that dephased population states are stationary under Az+Bz
result=test_close(result,'dephased singlet stationarity',Hz*Sd-Sd*Hz,0*Sd,1e-12,1e-12,...
                  'dephased deuteron-pair singlet population must commute with Az+Bz');
for n=1:numel(Td)
    result=test_close(result,['dephased triplet stationarity ' int2str(n)],Hz*Td{n}-Td{n}*Hz,0*Td{n},1e-12,1e-12,...
                      'dephased deuteron-pair triplet populations must commute with Az+Bz');
end
for n=1:numel(Qd)
    result=test_close(result,['dephased quintet stationarity ' int2str(n)],Hz*Qd{n}-Qd{n}*Hz,0*Qd{n},1e-12,1e-12,...
                      'dephased deuteron-pair quintet populations must commute with Az+Bz');
end

% Build a spherical-tensor system for captured stateinfo() output
sys_s.magnet=0;
sys_s.isotopes={'1H'};
inter_s.zeeman.scalar={0};
inter_s.temperature=300;
bas_s.formalism='sphten-liouv';
bas_s.approximation='none';
spin_s=test_spin_system(sys_s,inter_s,bas_s);
spin_s.sys.output=1;

% Capture and inspect the printed state composition report
printed=evalc('stateinfo(spin_s,state(spin_s,''Lz'',''1H''),1);');
result=test_true(result,'stateinfo norm output',contains(printed,'state vector 2-norm'),...
                 'stateinfo() must print the state norm through report()');
result=test_true(result,'stateinfo population output',contains(printed,'most populated basis states'),...
                 'stateinfo() must print the requested population table heading');
result=test_true(result,'stateinfo coefficient output',contains(printed,'+5.000e-01'),...
                 'stateinfo() must print the dominant state coefficient');

% Build a triplet-electron Hilbert-space system
sys_e.magnet=0;
sys_e.isotopes={'E3'};
inter_e.zeeman.scalar={0};
inter_e.temperature=300;
bas_e.formalism='zeeman-hilb';
bas_e.approximation='none';
spin_e=test_spin_system(sys_e,inter_e,bas_e);

% Project an isotropic zero-field population through arbitrary high-field axes
ZFS=[2 0.1 0.2; 0.1 1 0.3; 0.2 0.3 -3]*1e6;
Z=diag([27.8 28.1 28.4])*1e9;
rho_zf=zftrip(spin_e,ZFS,[1/3 1/3 1/3],Z,0.05,1);

% Check that an isotropic triplet population remains maximally mixed
result=test_close(result,'isotropic zftrip projection',rho_zf,speye(3)/3,1e-12,1e-12,...
                  'an isotropic zero-field triplet population is invariant under projection and dephasing');

end


