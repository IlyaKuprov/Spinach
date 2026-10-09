% Tests a one-spin liquid-state free induction decay. Syntax:
%
%                    result=test_liquid_single_spin_fid()
%
% Outputs:
%
%     result  - regression test result with explanatory messages
%
% The test simulates a zero-offset one-spin FID. With no Hamiltonian and no
% relaxation, transverse magnetisation is constant in time.
% The proton carrier is fixed at its original 14.1 T regression value,
% 3772062842.904 rad/s, independently of isotope-data revisions.
%
% ilya.kuprov@weizmann.ac.il

function result=test_liquid_single_spin_fid()

% Announce the test target
fprintf('TESTING: Single-spin liquid-state FID\n');

% State the NMR target of the test
result=new_test_result('kernel/liquid_single_spin_fid',...
                       'Single-spin liquid-state FID',...
                       'a zero-offset isolated spin has a constant free induction decay.');

% Build a one-spin Liouville-space system
sys.magnet=3772062842.904/spin('1H');
sys.isotopes={'1H'};
inter.zeeman.scalar={0};
bas.formalism='sphten-liouv';
bas.approximation={'none'};
spin_system=test_spin_system(sys,inter,bas);

% Set up a zero-offset acquisition
parameters.spins={'1H'};
parameters.rho0=state(spin_system,'L+','1H');
parameters.coil=state(spin_system,'L+','1H');
parameters.decouple={};
parameters.offset=0;
parameters.sweep=1000;
parameters.npoints=8;
parameters.zerofill=8;
parameters.axis_units='Hz';
parameters.invert_axis=0;

% Simulate the FID
fid=liquid(spin_system,@acquire,parameters,'nmr');

% The signal is constant and equal to its first point
result=test_close(result,'constant zero-offset FID',fid,fid(1)*ones(size(fid)),1e-12,1e-12,...
                  'without precession or relaxation the detected coherence is time-independent');

end


