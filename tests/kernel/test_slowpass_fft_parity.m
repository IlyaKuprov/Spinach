% Tests slowpass amplitude normalisation against time-domain FFT. Syntax:
%
%                    result=test_slowpass_fft_parity()
%
% Outputs:
%
%     result  - regression test result with explanatory messages
%
% The test compares frequency-domain acquisition from slowpass() with the
% same damped one-spin signal acquired in the time domain and processed by
% Matlab's unnormalised FFT.
%
% ilya.kuprov@weizmann.ac.il

function result=test_slowpass_fft_parity()

% Announce the test target
fprintf('TESTING: Slowpass FFT amplitude parity\n');

% State the slowpass normalisation target of the test
result=new_test_result('kernel/slowpass_fft_parity',...
                       'Slowpass FFT amplitude parity',...
                       'slowpass() must match the unnormalised amplitude convention of fft(acquire()).');

% Build a damped one-spin Liouville-space system
sys.magnet=14.1;
sys.isotopes={'1H'};
inter.zeeman.scalar={0};
inter.relaxation={'damp'};
inter.damp_rate=8.0;
inter.equilibrium='zero';
inter.rlx_keep='labframe';
inter.temperature=298;
bas.formalism='sphten-liouv';
bas.approximation='none';
spin_system=test_spin_system(sys,inter,bas);

% Get production generators and states
spin_system=assume(spin_system,'nmr');
H=hamiltonian(spin_system);
R=relaxation(spin_system);
K=kinetics(spin_system);
parameters.rho0=state(spin_system,'L+','1H');
parameters.coil=state(spin_system,'L+','1H');

% Acquire the signal in the time domain
parameters.decouple={};
parameters.sweep=4096;
parameters.npoints=4096;
fid=acquire(spin_system,parameters,H,R,K);
spectrum_fft=fftshift(fft(fid));

% Use the exact FFT bins as the slowpass frequency grid
frq_axis=ft_axis(0,parameters.sweep,parameters.npoints);
parameters.sweep=[frq_axis(1) frq_axis(end)];
spectrum_slow=slowpass(spin_system,parameters,H,R,K);

% Select the zero-frequency resonance
peak_idx=parameters.npoints/2+1;

% Check the FFT-compatible amplitude normalisation
result=test_close(result,'slowpass FFT amplitude',spectrum_slow(peak_idx),spectrum_fft(peak_idx),1e-6,2e-3,...
                  'frequency-domain acquisition should match the unnormalised FFT amplitude of the same damped FID');

% Check that single-point spectra are rejected before scaling
parameters_single=parameters;
parameters_single.sweep=[0 0];
parameters_single.npoints=1;
try
    slowpass(spin_system,parameters_single,H,R,K);
    single_rejected=false;
catch err
    single_rejected=contains(err.message,'integer greater than one');
end


result=test_true(result,'slowpass single-point grumbler',single_rejected,...
                 'single-point slowpass spectra cannot be normalised to the unnormalised FFT convention');

end

