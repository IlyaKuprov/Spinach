% Tests static powder and crystal contexts at one orientation. Syntax:
%
%                    result=test_ctx_powder_crystal()
%
% Outputs:
%
%     result  - regression test result with explanatory messages
%
% The test uses the single_crystal grid so that powder() and crystal()
% represent the same Euler orientation of the same anisotropic one-spin
% Hamiltonian.
%
% ilya.kuprov@weizmann.ac.il

function result=test_ctx_powder_crystal()

% Announce the test target
fprintf('TESTING: Powder and crystal single orientation\n');

% State the static-context target of the test
result=new_test_result('kernel/ctx_powder_crystal',...
                       'Powder and crystal single orientation',...
                       'powder() with the single_crystal grid must match crystal().');

% Build a one-spin anisotropic Liouville-space system
sys.magnet=14.1;
sys.isotopes={'1H'};
inter.zeeman.eigs={[-2 -2 4]};
inter.zeeman.euler={[0 0 0]};
bas.formalism='sphten-liouv';
bas.approximation='none';
bas.projections=+1;
spin_system=test_spin_system(sys,inter,bas);

% Set up a short static acquisition
parameters.spins={'1H'};
parameters.rho0=state(spin_system,'L+','1H');
parameters.coil=state(spin_system,'L+','1H');
parameters.decouple={};
parameters.offset=0;
parameters.sweep=2000;
parameters.npoints=4;
parameters.grid='single_crystal';
parameters.orientation=[0 0 0];
parameters.serial=true;
parameters.verbose=0;

% Run both production static contexts
fid_powder=powder(spin_system,@acquire,parameters,'nmr');
fid_crystal=crystal(spin_system,@acquire,parameters,'nmr');

% Check that one-point powder averaging is a crystal calculation
result=test_close(result,'single orientation powder',fid_powder,fid_crystal,1e-12,1e-12,...
                  'the single_crystal powder grid has unit weight at zero Euler angles');

end


