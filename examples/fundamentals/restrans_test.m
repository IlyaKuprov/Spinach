% Resonator transform test. Sends a square pulse into a simple
% resonator model and plots the time-domain response.
%
% ilya.kuprov@weizmann.ac.il
% u.rasulov@soton.ac.uk

function restrans_test()

% Get a pulse shape
pulse_in=randn(51,2);
pulse_in(1:3,:)=0;
pulse_in((end-2):end,:)=0;

% Get a figure going
figure(); scale_figure([2.0 1.0]);

% RLC response simulation, 1H @ 14.1 Tesla
subplot(1,2,1);
restrans(pulse_in(:,1), ... % X component of the pulse
         pulse_in(:,2), ... % Y component of the pulse
         1e-6,          ... % waveform slice duration
         2*pi*600e6,    ... % circuit resonance frequency
         50,            ... % circuit Q-factor
         'pwc',100);    ... % piecewise-constant model
ktitle('$^{1}$H @ 14.1 Tesla');

% RLC response simulation, 15N @ 14.1 Tesla
subplot(1,2,2);
restrans(pulse_in(:,1), ... % X component of the pulse
         pulse_in(:,2), ... % Y component of the pulse
         1e-6,          ... % waveform slice duration
         2*pi*60e6,     ... % circuit resonance frequency
         50,            ... % circuit Q-factor
         'pwc',100);    ... % piecewise-constant model
ktitle('$^{15}$N @ 14.1 Tesla');

end

