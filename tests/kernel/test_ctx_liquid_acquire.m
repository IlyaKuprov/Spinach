% Tests the liquid context against the direct acquire path. Syntax:
%
%                    result=test_ctx_liquid_acquire()
%
% Outputs:
%
%     result  - regression test result with explanatory messages
%
% The test runs a one-spin offset FID through liquid() and compares it
% with the same Hamiltonian, relaxation, and kinetics objects passed
% directly to acquire().
%
% ilya.kuprov@weizmann.ac.il

function result=test_ctx_liquid_acquire()

% Announce the test target
fprintf('TESTING: Liquid context acquire path\n');

% State the liquid-context target of the test
result=new_test_result('kernel/ctx_liquid_acquire',...
                       'Liquid context acquire path',...
                       'liquid() must pass the offset Liouvillian to acquire().');

% Build a one-spin Liouville-space system
sys.magnet=14.1;
sys.isotopes={'1H'};
inter.zeeman.scalar={0};
bas.formalism='sphten-liouv';
bas.approximation='none';
bas.projections=+1;
spin_system=test_spin_system(sys,inter,bas);

% Set up a short offset acquisition
parameters.spins={'1H'};
parameters.rho0=state(spin_system,'L+','1H');
parameters.coil=state(spin_system,'L+','1H');
parameters.decouple={};
parameters.offset=25;
parameters.sweep=1000;
parameters.npoints=4;
parameters.needs={};
parameters.rframes={};

% Run the production liquid context
fid_ctx=liquid(spin_system,@acquire,parameters,'nmr');

% Build the same direct acquire input path
spin_system=assume(spin_system,'nmr');
H=hamiltonian(spin_system);
R=relaxation(spin_system);
K=kinetics(spin_system);
H=frqoffset(spin_system,H,parameters);
fid_ref=acquire(spin_system,parameters,H,R,K);

% Check that the context passes through the direct dynamics
result=test_close(result,'liquid context FID',fid_ctx,fid_ref,1e-12,1e-12,...
                  'liquid() should reproduce direct acquire() for the same offset Liouvillian');

end


