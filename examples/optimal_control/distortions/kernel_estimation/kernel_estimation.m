% HiPER instrument filter function kernel estimation from 
% the quadrature components recorded by an antenna placed
% close to the sample location.
%
% Rob Hunter, Hassane el-Mkami, Graham Smith,
% Yujie Zhao, Shebha Anandhi Jegadeesan,
% Guinevere Mathies, Ilya Kuprov

function kernel_estimation()

% Read the experimental data
load('xix_on_resonance.mat','time_ns','real_part','imag_part');

% Plot the experimental data
kfigure(); scale_figure([1.5 2.0]);
subplot(3,1,1); plot(time_ns,[real_part imag_part]);
xlim tight; ylim padded; kgrid; kxlabel('time, ns');
kylabel('mV'); ktitle('HiPER instrument readout');

% Make ideal XiX waveform with ten 36 ns periods
xix_block=[ones(36,1); -ones(36,1)];
real_ideal=[zeros(175,1); xix_block; xix_block; xix_block; xix_block; xix_block; 
            xix_block; xix_block; xix_block; xix_block; xix_block; zeros(257,1)];
imag_ideal=zeros(1152,1);

% Plot the ideal waveform
subplot(3,1,2); plot(time_ns,[real_ideal imag_ideal]);
xlim tight; ylim padded; kgrid; kxlabel('time, ns');
kylabel('a.u.'); ktitle('ideal XiX waveform');

% Extract the kernel
x=real_ideal+1i*imag_ideal;
y=real_part+1i*imag_part;
h=kernelest(x,y,32,'tikh','causal',10);

% Compute and plot the convolution
y=conv(x,h); y=y(1:numel(x));
subplot(3,1,3); plot(time_ns,[real(y) imag(y)]);
xlim tight; ylim padded; kgrid; kxlabel('time, ns');
kylabel('mV'); ktitle('ideal + filter function');

% Plot the kernel
kfigure(); scale_figure([1.5 0.75]);
time_axis=0.5*(0:31)'; subplot(1,2,1); 
plot(time_axis,[real(h) imag(h)]);
xlim tight; ylim padded; kgrid;
ktitle('HiPER filter function kernel');
kxlabel('time, ns'); kylabel('mV');

% Plot the frequency response
h=abs(fftshift(fft(h,128)));
freq_axis=fft_freq_axis(32,0.5,128-32);
subplot(1,2,2); plot(freq_axis,h);
xlim tight; ylim padded; kgrid;
ktitle('HiPER frequency response');
kxlabel('frequency, GHz'); 
kylabel('abs. value, a.u.');

end

