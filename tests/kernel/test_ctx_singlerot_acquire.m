% Tests the single-rotor context with acquire(). Syntax:
%
%                    result=test_ctx_singlerot_acquire()
%
% Outputs:
%
%     result  - regression test result with explanatory messages
%
% The test runs a tiny anisotropic one-spin MAS calculation through
% singlerot() and checks the returned time-domain trace for basic
% physical and dimensional invariants.
%
% ilya.kuprov@weizmann.ac.il

function result=test_ctx_singlerot_acquire()

% State the single-rotor target of the test
result=new_test_result('kernel/ctx_singlerot_acquire',...
                       'Single-rotor acquire path',...
                       'singlerot() must project states into rotor space and run acquire().');

% Build a one-spin anisotropic Liouville-space system
sys.magnet=14.1;
sys.isotopes={'1H'};
inter.zeeman.eigs={[-2 -2 4]};
inter.zeeman.euler={[0 0 0]};
bas.formalism='sphten-liouv';
bas.approximation='none';
bas.projections=+1;
spin_system=test_spin_system(sys,inter,bas);

% Set up a tiny MAS acquisition
parameters.spins={'1H'};
parameters.rho0=state(spin_system,'L+','1H');
parameters.coil=state(spin_system,'L+','1H');
parameters.decouple={};
parameters.offset=0;
parameters.sweep=2000;
parameters.npoints=3;
parameters.rate=1000;
parameters.axis=[1 1 1];
parameters.max_rank=1;
parameters.grid='single_crystal';
parameters.serial=true;
parameters.verbose=0;

% Run the production single-rotor context
fid=singlerot(spin_system,@acquire,parameters,'nmr');

% Check the number of acquired points
result=test_close(result,'singlerot FID length',numel(fid),parameters.npoints,0,0,...
                  'acquire() should return one point for each requested time sample');

% Check that the zero-time signal survives state projection
fid_zero=parameters.coil'*parameters.rho0;
result=test_close(result,'singlerot zero-time signal',fid(1),fid_zero,1e-12,1e-12,...
                  'rotor-space projection must preserve the initial coil overlap');

% Check that the acquired trace is finite
result=test_true(result,'singlerot finite FID',all(isfinite(fid(:))),...
                 'short MAS propagation should not produce NaN or Inf values');

end


