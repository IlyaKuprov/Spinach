% Tests the double-rotor context with acquire(). Syntax:
%
%                    result=test_ctx_doublerot_acquire()
%
% Outputs:
%
%     result  - regression test result with explanatory messages
%
% The test runs a tiny anisotropic one-spin double-rotation calculation
% through doublerot() and checks the returned time-domain trace for basic
% physical and dimensional invariants.
%
% ilya.kuprov@weizmann.ac.il

function result=test_ctx_doublerot_acquire()

% State the double-rotor target of the test
result=new_test_result('kernel/ctx_doublerot_acquire',...
                       'Double-rotor acquire path',...
                       'doublerot() must project states into double-rotor space and run acquire().');

% Build a one-spin anisotropic Liouville-space system
sys.magnet=14.1;
sys.isotopes={'1H'};
inter.zeeman.eigs={[-2 -2 4]};
inter.zeeman.euler={[0 0 0]};
bas.formalism='sphten-liouv';
bas.approximation='none';
bas.projections=+1;
spin_system=test_spin_system(sys,inter,bas);

% Set up a tiny double-rotation acquisition
parameters.spins={'1H'};
parameters.rho0=state(spin_system,'L+','1H');
parameters.coil=state(spin_system,'L+','1H');
parameters.decouple={};
parameters.offset=0;
parameters.sweep=2000;
parameters.npoints=3;
parameters.rate_outer=800;
parameters.rate_inner=2400;
parameters.rank_outer=1;
parameters.rank_inner=1;
parameters.axis_outer=[sqrt(2/3) 0 sqrt(1/3)];
parameters.axis_inner=[sqrt(20-2*sqrt(30)) 0 sqrt(15+2*sqrt(30))];
parameters.grid='single_crystal';
parameters.serial=true;
parameters.verbose=0;

% Run the production double-rotor context
fid=doublerot(spin_system,@acquire,parameters,'nmr');

% Check the number of acquired points
result=test_close(result,'doublerot FID length',numel(fid),parameters.npoints,0,0,...
                  'acquire() should return one point for each requested time sample');

% Check that the zero-time signal survives double-rotor projection
fid_zero=parameters.coil'*parameters.rho0;
result=test_close(result,'doublerot zero-time signal',fid(1),fid_zero,1e-12,1e-12,...
                  'double-rotor projection must preserve the initial coil overlap');

% Check that the acquired trace is finite
result=test_true(result,'doublerot finite FID',all(isfinite(fid(:))),...
                 'short double-rotor propagation should not produce NaN or Inf values');

end


