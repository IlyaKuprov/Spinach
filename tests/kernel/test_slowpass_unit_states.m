% Slow passage spectrum with several substances: every substance has
% its own stationary unit state, and the frequency grid contains the
% exact zero frequency where the unit states used to make the linear
% solve singular. The slowpass spectrum must be finite on the whole
% grid and must agree with the FFT of the acquired FID at the peaks.
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=test_slowpass_unit_states.m>

function result=test_slowpass_unit_states()

% Test record
result=new_test_result('kernel/slowpass_unit_states',...
                       'Slowpass unit states per substance',...
                       'slowpass() must stay finite and FFT-consistent when several stationary unit states are present.');

% Two one-spin substances, no reactions, zero equilibrium
sys.magnet=14.1; sys.isotopes={'1H','1H'};
inter.chem.parts={1,2}; inter.chem.concs=[1 1];
inter.zeeman.scalar={0,0.05}; inter.coupling.scalar=cell(2);
inter.relaxation={'damp'}; inter.damp_rate=8.0;
inter.equilibrium='zero'; inter.rlx_keep='labframe';
inter.temperature=298;
bas.formalism='sphten-liouv'; bas.approximation={'none','none'};
spin_system=test_spin_system(sys,inter,bas);
spin_system=assume(spin_system,'nmr');
H=hamiltonian(spin_system); R=relaxation(spin_system);
K=kinetics(spin_system);

% Acquisition with a grid that contains zero frequency
parameters.rho0=state(spin_system,'L+','1H');
parameters.coil=state(spin_system,'L+','1H');
parameters.decouple={}; parameters.sweep=4096;
parameters.npoints=4096;
fid=acquire(spin_system,parameters,H,R,K);
spectrum_fft=fftshift(fft(fid));
frq_axis=ft_axis(0,parameters.sweep,parameters.npoints);
parameters.sweep=[frq_axis(1) frq_axis(end)];
spectrum_slow=slowpass(spin_system,parameters,H,R,K);

% Every grid point must be finite
result=test_true(result,'slowpass finite on a grid containing zero frequency',...
                 all(isfinite(spectrum_slow)),...
                 'stationary unit states must not make the frequency domain solve singular');

% Amplitude parity at the zero-frequency bin and at the second
% resonance; the FFT of a truncated FID carries a constant base-
% line of half the first point, which bounds the absolute error
zero_idx=parameters.npoints/2+1;
second_frq=1e-6*inter.zeeman.scalar{2}*spin_system.inter.basefrqs(2)/(2*pi);
[~,second_idx]=min(abs(frq_axis-second_frq)); baseline=abs(fid(1))/2;
result=test_close(result,'slowpass FFT amplitude at zero frequency',...
                  spectrum_slow(zero_idx),spectrum_fft(zero_idx),2*baseline,2e-3,...
                  'the unit state projection must leave the traceless spectrum unchanged');
result=test_close(result,'slowpass FFT amplitude at the second resonance',...
                  spectrum_slow(second_idx),spectrum_fft(second_idx),2*baseline,2e-3,...
                  'the second substance resonance must match the FFT of the same damped FID');

end

