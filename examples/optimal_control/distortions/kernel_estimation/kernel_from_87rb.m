% Transmitter and probe distortion kernel of a 400 MHz Bruker
% spectrometer fitted with a 4 mm Phoenix MAS probe, estimated
% from an oscilloscope recording of an eight-block XiX wave-
% form on the 87Rb channel.
%
% The wall clock record is heterodyned into the rotating frame,
% the carrier frequency is refined from the residual phase drift
% accumulated within each XiX half-block, the remaining constant
% phase is removed, and the FIR kernel is then obtained by line-
% ar least squares from the ideal and the measured waveforms.
%
% The heterodyne is an analytic signal demodulation, which is
% zero-phase: the rotating frame signal is not delayed with re-
% spect to the oscilloscope record.
%
% Calculation time: seconds
%
% Maxi Keitel, Ilya Kuprov

function kernel_from_87rb()

% Kernel sampling parameters
dt_kernel=1e-7; n_taps=40; t_kernel=dt_kernel*n_taps;

% XiX waveform parameters
n_blocks=8; t_block=10e-6;

% Start of the XiX waveform in the oscilloscope record
t_start=1.11e-5;

% Nominal carrier frequency of the 87Rb channel
carrier_freq=131.1468520e6;

% Fractional window inside each half-block, clear of the transitions
win_start=0.32; win_end=0.48;

% Read the oscilloscope record
load('87Rb_160W_XiX_10us_8repeats.mat','scope_trace','time_scope');
dt=time_scope(2)-time_scope(1);

% Plot the raw scope readout
kfigure(); scale_figure([2.0 1.5]);
subplot(2,2,1); plot(time_scope*1e6,scope_trace);
xlim padded; ylim padded; kgrid;
kxlabel('time, $\mu$s'); kylabel('scope reading, a.u.');
ktitle('wall clock');

% Segment holding the waveform and its ringdown
n_start=round(t_start/dt);
n_end=n_start+round((n_blocks*t_block+t_kernel)/dt);

% Time axis of the segment, relative to the start of the waveform
time_grid=dt*(0:(n_end-n_start))';

% Fitting windows inside each XiX half-block
n_low=round((win_start+0.5*(0:(2*n_blocks-1))')*t_block/dt)+1;
n_up=round((win_end+0.5*(0:(2*n_blocks-1))')*t_block/dt)+1;

% Heterodyne the record at the nominal carrier frequency
[real_bruk,imag_bruk]=heterodyne(dt,scope_trace,carrier_freq);
cplx_bruk=real_bruk(n_start:n_end)+1i*imag_bruk(n_start:n_end);

% Residual carrier offset from the phase slope in each half-block
freq_drift=zeros(2*n_blocks,1);
for n=1:(2*n_blocks)
    fit_win=n_low(n):n_up(n);
    phase_fit=polyfit(time_grid(fit_win),...
                      unwrap(angle(cplx_bruk(fit_win))),1);
    freq_drift(n)=phase_fit(1);
end

% Refine the carrier frequency, skipping the leading transient block
carrier_freq=carrier_freq-mean(freq_drift(2:end))/(2*pi);

% Heterodyne the record again at the refined carrier frequency
[real_bruk,imag_bruk]=heterodyne(dt,scope_trace,carrier_freq);
cplx_bruk=real_bruk(n_start:n_end)+1i*imag_bruk(n_start:n_end);

% Constant phase of the nominally positive XiX half-blocks
block_phase=zeros(n_blocks,1);
for n=1:n_blocks
    fit_win=n_low(2*n-1):n_up(2*n-1);
    block_phase(n)=angle(mean(cplx_bruk(fit_win)));
end

% Remove the constant phase and normalise the amplitude
cplx_bruk=exp(-1i*angle(mean(exp(1i*block_phase))))*cplx_bruk;
cplx_bruk=cplx_bruk/max(abs(cplx_bruk));

% Resample onto the kernel time grid
real_bruk=resample(real(cplx_bruk),round(1/dt_kernel),round(1/dt));
imag_bruk=resample(imag(cplx_bruk),round(1/dt_kernel),round(1/dt));
cplx_bruk=real_bruk+1i*imag_bruk;
time_grid=dt_kernel*(0:(numel(cplx_bruk)-1))';

% Ideal XiX waveform, silent after the last block
xix_ideal=(-1).^floor(time_grid/(0.5*t_block));
xix_ideal(time_grid>=n_blocks*t_block)=0;

% Plot the heterodyned signal and the ideal waveform
subplot(2,2,2); plot(time_grid*1e6,xix_ideal); hold on;
plot(time_grid*1e6,[real(cplx_bruk) imag(cplx_bruk)]); hold off;
xlim padded; ylim padded; kgrid;
kxlabel('time, $\mu$s'); kylabel('amplitude, a.u.');
klegend({'ideal','in-phase','quadrature'},'Location','Best');
ktitle('heterodyned signal');

% Extract the kernel and normalise the DC component
h=kernelest(xix_ideal,cplx_bruk,n_taps,'pinv');
h=h/sum(h);

% Plot the kernel
subplot(2,2,3); plot(dt_kernel*(0:(n_taps-1))*1e6,[real(h) imag(h)]);
xlim padded; ylim padded; kgrid;
kxlabel('time, $\mu$s'); kylabel('amplitude, a.u.');
klegend({'in-phase','quadrature'},'Location','Best');
ktitle('impulse response kernel');

% Plot the frequency response of the kernel
n_zf=256-n_taps; freq_axis=fft_freq_axis(n_taps,dt_kernel,n_zf);
subplot(2,2,4); plot(freq_axis*1e-6,abs(fftshift(fft(h,n_taps+n_zf))));
xlim padded; ylim padded; kgrid;
kxlabel('frequency, MHz'); kylabel('abs. value, a.u.');
ktitle('frequency response');

% Save the kernel
save('kernel_87rb.mat','h','dt_kernel');

end

