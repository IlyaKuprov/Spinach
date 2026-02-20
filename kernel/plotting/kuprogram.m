% Phase-sensitive spectrogram with an NMR-style phase referencing
% that evolves with the clock at each frequency line. Syntax:
%
%                    kuprogram(ux,uy,dt,win_len)
%
% Parameters:
%
%       ux      - pulse x-component on a uniform time grid
%
%       uy      - pulse y-component on a uniform time grid
%
%       dt      - time step in seconds
%
%       win_len - Gaussian window length in points
% 
%       gbc     - Gaussian blur parameter for chroma
%
% Output:
%
%       a spectrogram that uses hue to encode phase
%       and brightness to encode amplitude
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=kuprogram.m>

function kuprogram(ux,uy,dt,win_len,gbc)

% Check consistency
grumble(ux,uy,dt,win_len,gbc);

% Inputs to a complex column
ux=ux(:); uy=uy(:); pulse=ux+1i*uy;

% Sampling frequency
samp_freq=1/dt;

% Amp-normalised Gaussian window
sam_idx=(0:win_len-1).';
cen_samp=(win_len-1)/2;
sig_samp=(win_len-1)/6;
win_fun=exp(-0.5*((sam_idx-cen_samp)/sig_samp).^2);
win_fun=win_fun/sum(win_fun);

% Hop size and the FFT size
hop=max(1,round(win_len/4));
nfft=2^nextpow2(win_len);

% Determine the number of samples in the pulse
n_samp=numel(pulse);

% Pad the pulse to fit an integer number of frames
n_frm=1+ceil((n_samp-win_len)/hop);
pad_len=(n_frm-1)*hop+win_len;
if(pad_len>n_samp)
    pulse=[pulse; zeros(pad_len-n_samp,1)];
end

% Windowed FFT for each frame
spec=complex(zeros(nfft,n_frm));
t_start=zeros(1,n_frm);
for k=1:n_frm
    idx0=(k-1)*hop;
    seg=pulse(idx0+1:idx0+win_len).*win_fun;
    spec(:,k)=fft(seg,nfft);
    t_start(k)=idx0*dt;
end

% Frame time axis at window centres
t_ctr=t_start+((win_len-1)/2)*dt;

% Centred frequency axis
if(mod(nfft,2)==0)
    freq_axis=(-nfft/2:nfft/2-1).'*samp_freq/nfft;
else
    freq_axis=(-(nfft-1)/2:(nfft-1)/2).'*samp_freq/nfft;
end

% Apply FFT shift to centre the spectrum
spec=fftshift(spec,1);

% Apply NMR-style phase referencing with s0=0
phs_cor=exp(-1i*2*pi*(freq_axis*t_start));
spec=spec.*phs_cor;

% Compute amplitude and phase maps
amp_map=abs(spec); phs_map=angle(spec);

% Map amplitude into brightness on a linear scale
amp_max=max(amp_map(:));
if(amp_max==0)
    amp_max=1;
end
val_map=amp_map/amp_max;
val_map=min(max(val_map,0),1);

% Map phase into hue with full saturation
hue_map=mod(phs_map,2*pi)/(2*pi);
sat_map=ones(size(hue_map));

% Convert HSV image into RGB
hsv_img=cat(3,hue_map,sat_map,val_map);
rgb_img=hsv2rgb(hsv_img);

% Blur the chroma channel
rgb_img_gb=imgaussfilt(rgb_img,gbc);
rgb_img_gb=rgb_img_gb/max(rgb_img_gb,[],'all');
ycbcr_img_gb=rgb2ycbcr(rgb_img_gb);
ycbcr_img=rgb2ycbcr(rgb_img);
ycbcr_img(:,:,2)=ycbcr_img_gb(:,:,2);
ycbcr_img(:,:,3)=ycbcr_img_gb(:,:,3);
rgb_img=ycbcr2rgb(ycbcr_img);

% Create the figure and display the phase-sensitive spectrogram
kfigure('Color','w');
image('XData',[t_ctr(1)     t_ctr(end)],...
      'YData',[freq_axis(1) freq_axis(end)],...
      'CData',rgb_img);
set(gca,'YDir','normal');

% Add axis labels and limits
kxlabel('time (s)');
kylabel('frequency (Hz)');
xlim([0,(n_samp-1)*dt]);

end

% Consistency enforcement
function grumble(ux,uy,dt,win_len,gbc)
if (~isnumeric(ux))||(~isreal(ux))||...
   (~isvector(ux))||isempty(ux)
    error('ux must be a non-empty real numeric vector');
end
if (~isnumeric(uy))||(~isreal(uy))||...
   (~isvector(uy))||isempty(uy)
    error('uy must be a non-empty real numeric vector');
end
if (numel(ux)~=numel(uy))
    error('ux and uy must have the same number of elements');
end
if (~isnumeric(dt))||(~isreal(dt))||(~isscalar(dt))||...
   (~isfinite(dt))||(dt<=0)
    error('dt must be a positive real scalar');
end
if (~isnumeric(win_len))||(~isreal(win_len))||...
   (~isscalar(win_len))||(~isfinite(win_len))||(win_len<=0)
    error('win_len must be a positive real scalar');
end
if (win_len~=round(win_len))
    error('win_len must be an integer');
end
if (win_len<7)
    error('win_len must be at least 7');
end
if (win_len>numel(ux))
    error('win_len must not exceed the number of samples in ux and uy');
end
if (~isnumeric(gbc))||(~isreal(gbc))||(~isscalar(gbc))||...
   (~isfinite(gbc))||(gbc<0)
    error('gbc must be a non-negative real scalar');
end
end

% "Email! The way all modern tragedies begin."
%
% Jessi Jezewska Stevens

