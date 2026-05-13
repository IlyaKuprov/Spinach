% Tests the grid-free Fokker-Planck context with acquire(). Syntax:
%
%                    result=test_ctx_gridfree_acquire()
%
% Outputs:
%
%     result  - regression test result with explanatory messages
%
% The test runs a tiny anisotropic one-spin MAS calculation through
% gridfree() and checks the returned time-domain trace for basic physical
% and dimensional invariants.
%
% ilya.kuprov@weizmann.ac.il

function result=test_ctx_gridfree_acquire()

% State the grid-free target of the test
result=new_test_result('kernel/ctx_gridfree_acquire',...
                       'Grid-free acquire path',...
                       'gridfree() must project states into SLE space and run acquire().');

% Build a one-spin anisotropic Liouville-space system
sys.magnet=14.1;
sys.isotopes={'1H'};
inter.zeeman.eigs={[-2 -2 4]};
inter.zeeman.euler={[0 0 0]};
bas.formalism='sphten-liouv';
bas.approximation='none';
bas.projections=+1;
spin_system=test_spin_system(sys,inter,bas);

% Set up a tiny grid-free acquisition
parameters.spins={'1H'};
parameters.rho0=state(spin_system,'L+','1H');
parameters.coil=state(spin_system,'L+','1H');
parameters.decouple={};
parameters.offset=0;
parameters.sweep=2000;
parameters.npoints=3;
parameters.rate=1000;
parameters.axis=[1 1 1];
parameters.max_rank=2;
parameters.verbose=0;

% Run the production grid-free context
fid=gridfree(spin_system,@acquire,parameters,'nmr');

% Check the number of acquired points
result=test_close(result,'gridfree FID length',numel(fid),parameters.npoints,0,0,...
                  'acquire() should return one point for each requested time sample');

% Check that the zero-time signal survives SLE projection
fid_zero=parameters.coil'*parameters.rho0;
result=test_close(result,'gridfree zero-time signal',fid(1),fid_zero,1e-12,1e-12,...
                  'SLE-space projection must preserve the initial coil overlap');

% Check that the acquired trace is finite
result=test_true(result,'gridfree finite FID',all(isfinite(fid(:))),...
                 'short grid-free propagation should not produce NaN or Inf values');

end


