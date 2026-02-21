% Response function extraction from the transmission profile
% of the HiPER instrument. We are sending a linear chirp via
% the AWG from 93.3 GHz to 94.7 GHz over 1400 ns and record-
% ing the power response (so a square root must be taken to
% obtain the amplitude). The two spikes show when the chirp
% starts and finishes.
%
% Rob Hunter, Hassane el-Mkami, Graham Smith,
% Yujie Zhao, Shebha Anandhi Jegadeesan,
% Guinevere Mathies, Ilya Kuprov

function kernel_from_trans()

% Load power data and convert to amplitude
load('power_at_eik.mat','freq_axis_ghz',...
                        'power_at_eik');
power_at_eik(power_at_eik<0)=0;
amp=sqrt(power_at_eik); amp=amp/max(amp);

% Plot as received
kfigure(); scale_figure([0.75 1.75]);
subplot(3,1,1); plot(freq_axis_ghz,amp);
kxlabel('frequency, GHz'); kgrid;
kylabel('filter ampl.');
ktitle('transmission spectrum');
xlim tight; ylim padded;

% Get apodisation weights
leave_intact=(freq_axis_ghz>93.4)&...
             (freq_axis_ghz<94.6);
np_intact=nnz(leave_intact);
fade_in=sin(linspace(0,pi/2,nnz(~leave_intact)/2)).^4;
weights=[fade_in ones(1,np_intact) fliplr(fade_in)];

% Apply weights and resample
amp=resample(amp.*weights',1,1000);
freq_axis_ghz=linspace(freq_axis_ghz(1),...
                       freq_axis_ghz(end),numel(amp))';
subplot(3,1,2); plot(freq_axis_ghz,amp); 
hold on; plot(freq_axis_ghz,...
              weights(1:1000:end),'r-');
kxlabel('frequency, GHz'); kgrid;
kylabel('filter ampl.');
ktitle('resampling and apodisation');
xlim tight; ylim padded;
klegend({'spectrum','window func.'},...
        'Location','South');

% Get the convolution kernel with a 0.5 ns time step
h=ifft(ifftshift([zeros(size(amp)); zeros(size(amp)); amp;
                  zeros(size(amp)); zeros(size(amp))]));
df=freq_axis_ghz(2)-freq_axis_ghz(1); zf=2*numel(amp);
[~,t]=ifft_time_axis(numel(amp),df,zf);
h_old=h(t<=16); t_old=t(t<=16); t=(0:31)'/2;
h=interp1(t_old,h_old,t,'spline');
save('hiper_kernel.mat','h');
subplot(3,1,3); plot(t,[real(h) imag(h)]);
kxlabel('time, ns'); kgrid; xlim tight; ylim padded;
kylabel('filter ampl.'); klegend({'real','imag'});
ktitle('filter kernel at 94.0 GHz offset');

end

