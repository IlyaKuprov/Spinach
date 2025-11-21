% An illustration of the effect of the resonator response function
% on a typical composite pulse in NMR spectroscopy. The argument
% may be set to 'previous' (default, corresponds to piecewise con-
% stant input waveform) or any of the options ('linear', 'cubic',
% etc.) supported by interp1() function.
%
% Calculation time: seconds.
%
% ilya.kuprov@weizmann.ac.il
% u.rasulov@soton.ac.uk

function rlc_response_1(interp_type)

% Settings (14N NMR @ 14.09 Tesla)
omega=14.09*spin('14N'); Q=80; 
n_slices=20; tmax=50e-6;

% Default to piecewise-constant
if ~exist('interp_type','var')
    interp_type='previous';
end

% Decide the time grid (4 x Nyquist)
dt=pi/(4*omega); npts=tmax/dt+1;
time_grid=linspace(0,tmax,npts);

% Random amplitude component
amp_part=rand(1,n_slices+1);
amp_part(1:2)=0; amp_part((end-2):end)=0;
amp_part=interp1(linspace(0,tmax,n_slices+1),...
                 amp_part,time_grid,interp_type);

% Random phase component
phi_part=2*pi*rand(1,n_slices+1);
phi_part(1:2)=0; phi_part((end-2):end)=0;
phi_part=interp1(linspace(0,tmax,n_slices+1),...
                 phi_part,time_grid,interp_type);

% Put the pulse together
inp_signal=amp_part.*cos(omega*time_grid+phi_part);

% Heterodyne out the carrier frequency
inp_real=amp_part.*cos(phi_part);
inp_imag=amp_part.*sin(phi_part);

% Plot input signal components
kfigure(); scale_figure([1.75 1.2]);
subplot(2,2,1); plot(1e6*time_grid,inp_signal);
kxlabel('time, $\mu$s'); kylabel('voltage, a.u.');
ktitle('input, wall clock'); 
kgrid; axis tight; ylim([-1.1 1.1]);
subplot(2,2,2); plot(1e6*time_grid,[inp_real; inp_imag]);
kxlabel('time, $\mu$s'); kylabel('voltage, a.u.');
ktitle('input, heterodyne at $\omega_{0}$');
kgrid; axis tight; ylim([-1.1 1.1]);
legend({'real','imag'},'Location','NorthEast');

% Build the RLC bandpass response kernel
sys=tf([1/(omega*Q) 0],[1/(omega^2) 1/(omega*Q) 1]);

% Apply the RLC bandpass response kernel
[out_signal,time_grid]=lsim(sys,inp_signal,time_grid);

% Heterodyne out the carrier frequency
out_real=+2*lowpass(out_signal.*cos(omega*time_grid),0.1);
out_imag=-2*lowpass(out_signal.*sin(omega*time_grid),0.1);

% Plot output signal components
subplot(2,2,3); plot(1e6*time_grid,out_signal);
kxlabel('time, $\mu$s'); kylabel('voltage, a.u.');
ktitle('output, wall clock'); 
kgrid; axis tight; ylim([-1.1 1.1]);
subplot(2,2,4); plot(1e6*time_grid,[out_real out_imag]);
kxlabel('time, $\mu$s'); kylabel('voltage, a.u.');
ktitle('output, heterodyne at $\omega_{0}$'); 
kgrid; axis tight; ylim([-1.1 1.1]);
legend({'real','imag'},'Location','NorthEast');

end

