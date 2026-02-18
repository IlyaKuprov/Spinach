% HiPER instrument filter function kernel application to 
% a complicated shaped pulse and a comparison with expe-
% rimental measurement at the instrument.
%
% Rob Hunter, Hassane el-Mkami, Graham Smith,
% Yujie Zhao, Shebha Anandhi Jegadeesan,
% Guinevere Mathies, Ilya Kuprov

function kernel_application()

% Read the input pulse
load('shaped_pulse_inp.mat','time_ns',...
     'real_part','imag_part');

% Plot the input pulse
kfigure(); scale_figure([1.0 1.5]);
subplot(3,1,1); plot(time_ns,[real_part imag_part]);
xlim tight; ylim padded; kgrid; kxlabel('time, ns');
kylabel('a.u.'); ktitle('input pulse');

% Load HiPER kernel
load('hiper_kernel.mat','h');

% Compute and plot the convolution
x=real_part+1i*imag_part;
y=conv(x,h); y=y(1:numel(x));
subplot(3,1,2); plot(time_ns,[real(y) imag(y)]);
xlim tight; ylim padded; kgrid; kxlabel('time, ns');
kylabel('a.u.'); ktitle('ideal $+$ filter function');

% Read the measured pulse
load('shaped_pulse_out.mat','time_ns',...
     'real_part','imag_part');
real_part=real_part(174:934);
imag_part=imag_part(174:934);
time_ns=time_ns(174:934);
time_ns=time_ns-time_ns(1);

% Plot the measured pulse
subplot(3,1,3); plot(time_ns,[real_part imag_part]);
xlim tight; ylim padded; kgrid; kxlabel('time, ns');
kylabel('a.u.'); ktitle('measured pulse');

end