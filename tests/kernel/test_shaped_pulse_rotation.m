% Tests a one-slice Cartesian shaped pulse. Syntax:
%
%                    result=test_shaped_pulse_rotation()
%
% Outputs:
%
%     result  - regression test result with explanatory messages
%
% The test applies one rectangular X pulse slice with amplitude 1 rad/s and
% duration pi seconds; the net flip angle is pi, so Lz must invert.
%
% ilya.kuprov@weizmann.ac.il

function result=test_shaped_pulse_rotation()

% State the pulse target of the test
result=new_test_result('kernel/shaped_pulse_rotation',...
                       'Piecewise Cartesian pulse rotation',...
                       'a rectangular Cartesian shaped pulse must reproduce the hard-pulse limit.');

% Build a one-proton Hilbert-space spin system
sys.magnet=0;
sys.isotopes={'1H'};
inter.zeeman.scalar={0};
bas.formalism='zeeman-hilb';
bas.approximation='none';
spin_system=test_spin_system(sys,inter,bas);

% Define drift, controls, and one pi pulse slice
Lx=operator(spin_system,'Lx',1);
Ly=operator(spin_system,'Ly',1);
Lz=state(spin_system,'Lz',1);
drift=0*Lx;
controls={Lx,Ly};
amplitudes={1,0};
slice_durs=pi;

% Apply the shaped pulse
rho_obs=shaped_pulse_xy(spin_system,drift,controls,amplitudes,slice_durs,Lz,'expm-pwc');

% Check the physical rotation result
result=test_close(result,'rectangular X pi pulse',rho_obs,-Lz,1e-14,1e-14,...
                  'amplitude times duration is the flip angle in radians');

end

