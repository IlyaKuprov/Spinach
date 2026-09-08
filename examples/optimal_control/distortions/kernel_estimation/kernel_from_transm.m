% Response function extraction from the transmission profile
% of the HiPER instrument. We are sending a linear chirp via
% the AWG from 93.3 GHz to 94.7 GHz over 1400 ns and record-
% ing the power response (so a square root must be taken to
% obtain the amplitude). The two spikes show when the chirp
% starts and finishes.
%
% The measurement has no phase information; the kernel is
% therefore taken to be the minimum-phase kernel with the
% measured amplitude spectrum, obtained by the real cepstrum
% (Kramers-Kronig) construction. The zero-phase alternative
% would be symmetric in time and thus non-causal.
%
% Rob Hunter, Hassane el-Mkami, Graham Smith,
% Yujie Zhao, Shebha Anandhi Jegadeesan,
% Guinevere Mathies, Ilya Kuprov

function kernel_from_transm()

% Load power data and convert to amplitude
load('power_at_eik.mat','freq_axis_ghz',...
                        'power_at_eik');
power_at_eik(power_at_eik<0)=0;
amp=sqrt(power_at_eik); amp=amp/max(amp);

% Plot as received
kfigure(); scale_figure([0.75 1.00]);
subplot(2,1,1); plot(freq_axis_ghz,amp,'Color',[0.85 0.85 0.85]);
hold on; kxlabel('frequency, GHz'); kgrid; xlim tight;
kylabel('filter ampl.'); ylim padded;
ktitle('transmission spectrum');

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
plot(freq_axis_ghz,weights(1:1000:end),'r-');
plot(freq_axis_ghz,amp,'k-'); 
klegend({'spectrum raw','window function',...
         'spectrum filtered'},'Location','Best');

% Zero-fill the amplitude spectrum and put 94.0 GHz at zero frequency
amp_spec=ifftshift([zeros(size(amp)); zeros(size(amp)); amp;
                    zeros(size(amp)); zeros(size(amp))]);

% Real cepstrum of the log-amplitude, floored at -60 dB
cep=ifft(log(max(amp_spec,1e-3)));

% Minimum-phase kernel by causal cepstrum folding, even-length grid
npts=numel(amp_spec); cep_win=[1; 2*ones(npts/2-1,1); 1; zeros(npts/2-1,1)];
h=ifft(exp(fft(cep_win.*cep)));

% Resample to a 0.5 ns time step and rescale to keep the DC gain
df=freq_axis_ghz(2)-freq_axis_ghz(1); zf=2*numel(amp);
[~,t,dt]=ifft_time_axis(numel(amp),df,zf);
h_old=h(t<=16)*(0.5/dt); t_old=t(t<=16); t=(0:31)'/2;
h=interp1(t_old,h_old,t,'spline');
save('hiper_kernel_trans.mat','h');
subplot(2,1,2); plot(t,[real(h) imag(h)]);
kxlabel('time, ns'); kgrid; xlim tight; ylim padded;
kylabel('filter ampl.'); klegend({'in-phase','quadrature'});
ktitle('filter kernel at 94.0 GHz offset');

end

