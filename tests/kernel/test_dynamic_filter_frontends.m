% Tests dynamic state-filter front-end kernels on compact spin systems. Syntax:
%
%                    result=test_dynamic_filter_frontends()
%
% Outputs:
%
%     result  - regression test result with explanatory messages
%
% The test exercises coherence(), correlation(), decouple(), homospoil(),
% and spinlock() on small spherical-tensor Liouville-space systems with
% analytically known surviving state components.
%
% ilya.kuprov@weizmann.ac.il

function result=test_dynamic_filter_frontends()

% State the dynamic filter target of the test
result=new_test_result('kernel/dynamic_filter_frontends',...
                       'Dynamic state-filter front ends',...
                       'state-selection front ends must keep exactly the intended basis components.');

% Check analytical coherence and correlation filters
result=local_test_coherence_correlation(result);

% Check analytical decoupling of a coupled spin system
result=local_test_decoupling(result);

% Check homospoil zero-quantum and longitudinal filters
result=local_test_homospoil(result);

% Check analytical spin-lock projection along both transverse axes
result=local_test_spinlock(result);

end


function result=local_test_coherence_correlation(result)

% Build a heteronuclear spherical-tensor Liouville-space system
sys.magnet=14.1;
sys.isotopes={'1H','13C'};
inter.zeeman.scalar={1.0,2.0};
inter.coupling.scalar=cell(2);
inter.coupling.scalar{1,2}=10.0;
inter.coupling.scalar{2,2}=0.0;
bas.formalism='sphten-liouv';
bas.approximation='none';
spin_system=test_spin_system(sys,inter,bas);

% Construct physical single- and two-spin state components
rho_hz=state(spin_system,'Lz','1H');
rho_cz=state(spin_system,'Lz','13C');
rho_hp_cz=state(spin_system,{'L+','Lz'},{1,2});
rho_hm_cz=state(spin_system,{'L-','Lz'},{1,2});
rho_hz_cz=state(spin_system,{'Lz','Lz'},{1,2});

% Select only the proton plus-one coherence from a mixed stack
rho_stack=[rho_hp_cz+rho_hm_cz 2*rho_hp_cz+3*rho_hm_cz];
rho_keep=[rho_hp_cz 2*rho_hp_cz];
rho_obs=coherence(spin_system,rho_stack,{{'1H',1}});
result=test_close(result,'coherence isotope selection',rho_obs,rho_keep,1e-14,1e-14,...
                  'coherence() must retain only the requested proton coherence order');

% Select only one-spin correlations from all spins
rho_mix=rho_hz+2*rho_cz+3*rho_hz_cz;
rho_obs=correlation(spin_system,rho_mix,1,'all');
rho_keep=rho_hz+2*rho_cz;
result=test_close(result,'correlation all-spin selection',rho_obs,rho_keep,1e-14,1e-14,...
                  'correlation() must retain one-spin terms and reject two-spin order');

% Select only carbon-local terms using the numeric spin-list path
rho_obs=correlation(spin_system,rho_mix,1,2);
rho_keep=2*rho_cz+3*rho_hz_cz;
result=test_close(result,'correlation numeric spin selection',rho_obs,rho_keep,1e-14,1e-14,...
                  'correlation() must honour numeric spin lists when counting correlations');

end


function result=local_test_decoupling(result)

% Build a coupled heteronuclear system and request the NMR Hamiltonian
spin_system=local_heteronuclear_system();
spin_system=assume(spin_system,'nmr');
H=hamiltonian(spin_system);

% Build a state with carbon-local and proton-carbon terms present
rho_h=state(spin_system,'Lz','1H');
rho_c=state(spin_system,'Lz','13C');
rho_hc=state(spin_system,{'Lz','Lz'},{1,2});
rho=rho_h+2*rho_c+3*rho_hc;

% Decouple carbon by isotope name
[H_obs,rho_obs]=decouple(spin_system,H,rho,{'13C'});
result=test_close(result,'decoupled state by name',rho_obs,rho_h,1e-14,1e-14,...
                  'decouple() must remove every state component involving the named spin');

% Check that rows and columns touching carbon have been removed
carbon_mask=(spin_system.bas.basis(:,2)~=0);
zero_rows=norm(H_obs(carbon_mask,:),1);
zero_cols=norm(H_obs(:,carbon_mask),1);
result=test_close(result,'decoupled Liouvillian rows',zero_rows+zero_cols,0,1e-14,1e-14,...
                  'decouple() must zero Liouvillian rows and columns involving the decoupled spin');

% Decouple carbon by numeric spin index
[~,rho_obs]=decouple(spin_system,H,rho,2);
result=test_close(result,'decoupled state by index',rho_obs,rho_h,1e-14,1e-14,...
                  'decouple() must accept numeric spin indices as well as isotope names');

end


function result=local_test_homospoil(result)

% Build a homonuclear two-spin system with genuine zero-quantum states
sys.magnet=14.1;
sys.isotopes={'1H','1H'};
inter.zeeman.scalar={0.0,0.0};
inter.coupling.scalar=cell(2);
inter.coupling.scalar{1,2}=0.0;
inter.coupling.scalar{2,2}=0.0;
bas.formalism='sphten-liouv';
bas.approximation='none';
spin_system=test_spin_system(sys,inter,bas);

% Mix longitudinal, zero-quantum, and double-quantum components
rho_long=state(spin_system,{'Lz','Lz'},{1,2});
rho_zq=state(spin_system,{'L+','L-'},{1,2});
rho_dq=state(spin_system,{'L+','L+'},{1,2});
rho_mix=rho_long+2*rho_zq+3*rho_dq;

% Keep zero-quantum coherences under the experimental approximation
rho_obs=homospoil(spin_system,rho_mix,'keep');
rho_keep=rho_long+2*rho_zq;
result=test_close(result,'homospoil zero-quantum keep',rho_obs,rho_keep,1e-14,1e-14,...
                  'homospoil() keep mode must retain zero-quantum coherences and longitudinal terms');

% Destroy zero-quantum coherences for the strict longitudinal projection
rho_obs=homospoil(spin_system,rho_mix,'destroy');
result=test_close(result,'homospoil zero-quantum destroy',rho_obs,rho_long,1e-14,1e-14,...
                  'homospoil() destroy mode must retain only longitudinal spherical tensors');

end


function result=local_test_spinlock(result)

% Build a one-spin system for analytical transverse locking
sys.magnet=14.1;
sys.isotopes={'1H'};
inter.zeeman.scalar={0.0};
bas.formalism='sphten-liouv';
bas.approximation='none';
spin_system=test_spin_system(sys,inter,bas);
spin_system=assume(spin_system,'nmr');

% Prepare Cartesian generators and a mixed transverse-longitudinal state
Lx=operator(spin_system,'Lx','1H');
Ly=operator(spin_system,'Ly','1H');
rho_x=state(spin_system,'Lx','1H');
rho_y=state(spin_system,'Ly','1H');
rho_z=state(spin_system,'Lz','1H');
rho_mix=rho_x+2*rho_y+3*rho_z;

% Lock along X and compare against the retained X magnetisation
rho_obs=spinlock(spin_system,Lx,Ly,rho_mix,'X');
result=test_close(result,'spinlock X projection',rho_obs,rho_x,1e-13,1e-13,...
                  'spinlock() along X must destroy Y and Z components in the locked frame');

% Lock along Y and compare against the retained Y magnetisation
rho_obs=spinlock(spin_system,Lx,Ly,rho_mix,'Y');
result=test_close(result,'spinlock Y projection',rho_obs,2*rho_y,1e-13,1e-13,...
                  'spinlock() along Y must destroy X and Z components in the locked frame');

end


function spin_system=local_heteronuclear_system()

% Build the common heteronuclear spherical-tensor test system
sys.magnet=14.1;
sys.isotopes={'1H','13C'};
inter.zeeman.scalar={1.0,2.0};
inter.coupling.scalar=cell(2);
inter.coupling.scalar{1,2}=10.0;
inter.coupling.scalar{2,2}=0.0;
bas.formalism='sphten-liouv';
bas.approximation='none';
spin_system=test_spin_system(sys,inter,bas);

end


