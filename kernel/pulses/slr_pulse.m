% Shinnar-Le Roux linear-phase selective excitation pulse. Syntax:
%
% [Cx,Cy,durs,amps,phis]=slr_pulse(npts,dur,tbw,flip_angle,pass_rip,stop_rip)
%
% Parameters:
%
%       npts - even number of piecewise-constant pulse slices
%
%        dur - total pulse duration, seconds
%
%        tbw - time-bandwidth product, defined as the pulse
%              duration times the nominal full passband width
%
% flip_angle - on-resonance flip angle between zero and pi/2,
%              radians
%
%   pass_rip - 90-degree excitation passband ripple target used
%              in the prototype design, dimensionless
%
%   stop_rip - 90-degree excitation stopband ripple target used
%              in the prototype design, dimensionless
%
% Outputs:
%
%         Cx - X control amplitudes, rad/s, 1 x npts row vector
%
%         Cy - Y control amplitudes, rad/s, 1 x npts row vector
%
%       durs - pulse slice durations, seconds, 1 x npts row vector
%
%       amps - RF amplitudes, rad/s, 1 x npts row vector
%
%       phis - RF phases, radians, 1 x npts row vector
%
% The beta polynomial is obtained by continuous weighted least squares
% in a linear-phase cosine basis. The complementary minimum-phase alpha
% polynomial and the RF waveform are then obtained by the inverse SLR
% transform. The ripple arguments enter the excitation-pulse transform
% and the transition-width estimate of Pauly et al.; they are design
% targets rather than guaranteed minimax error bounds.
% For flip angles below pi/2, they do not specify angle-independent
% magnetisation error bounds.
%
% The output controls are calibrated for Spinach propagation under
% exp(-1i*H*t) and may be passed directly to shaped_pulse_xy().
%
% J. Pauly, P. Le Roux, D. Nishimura, and A. Macovski,
% IEEE Transactions on Medical Imaging 10(1), 53-65 (1991),
% https://doi.org/10.1109/42.75611
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=slr_pulse.m>

function [Cx,Cy,durs,amps,phis]=slr_pulse(npts,dur,tbw,flip_angle,pass_rip,stop_rip)

% Check consistency
grumble(npts,dur,tbw,flip_angle,pass_rip,stop_rip);

% Convert 90-degree excitation-profile targets into beta targets
beta_pass=sqrt(pass_rip/2);
beta_stop=stop_rip/sqrt(2);
beta_scale=sin(flip_angle/2);

% Estimate the transition width using the Pauly relation
log_pass=log10(beta_pass);
log_stop=log10(beta_stop);
trans_measure=(5.309e-3*log_pass^2+7.114e-2*log_pass-4.761e-1)*log_stop+...
              (-2.66e-3*log_pass^2-5.941e-1*log_pass-4.278e-1);
pass_edge=(tbw-trans_measure)/npts;
stop_edge=(tbw+trans_measure)/npts;

% Reject filter specifications that do not fit below Nyquist
if (pass_edge<=0)||(pass_edge>=stop_edge)||(stop_edge>=1)
    error('the specified sampling, bandwidth, and ripples do not produce a valid transition band.');
end

% Build the exact weighted least-squares normal equations
half_n=npts/2;
orders=(0:(half_n-1))+0.5;
diff_orders=orders'-orders;
sum_orders=orders'+orders;
pass_lim=pi*pass_edge;
stop_lim=pi*stop_edge;
zero_mask=(diff_orders==0);
pass_diff=sin(diff_orders*pass_lim)./diff_orders;
pass_diff(zero_mask)=pass_lim;
pass_sum=sin(sum_orders*pass_lim)./sum_orders;
stop_diff=(sin(diff_orders*pi)-sin(diff_orders*stop_lim))./diff_orders;
stop_diff(zero_mask)=pi-stop_lim;
stop_sum=(sin(sum_orders*pi)-sin(sum_orders*stop_lim))./sum_orders;
stop_weight=beta_pass/beta_stop;
gram=0.5*(pass_diff+pass_sum)+0.5*stop_weight*(stop_diff+stop_sum);
moment=sin(orders'*pass_lim)./orders';

% Solve the positive-definite least-squares system
gram=(gram+gram')/2;
[chol_factor,chol_flag]=chol(gram,'lower');
if chol_flag~=0
    error('the SLR least-squares filter design is numerically singular.');
end
cos_coeff=chol_factor'\(chol_factor\moment);
left_half=flipud(cos_coeff)/2;
beta=[left_half;flipud(left_half)]';

% Set the beta-polynomial centre response to the requested flip
centre_resp=sum(beta);
if abs(centre_resp)<=eps(norm(beta,1))
    error('the SLR least-squares filter has zero centre response.');
end
beta=(beta_scale/centre_resp)*beta;

% Sample the beta response on an oversampled Fourier grid
fft_len=2^(nextpow2(npts)+4);
beta_pad=zeros(1,fft_len);
beta_pad(1:npts)=beta;
beta_resp=fft(beta_pad);
resp_peak=max(abs(beta_resp));

% Require a strictly positive complementary spectrum
if resp_peak>=1
    error('the requested SLR beta response reaches or exceeds unit magnitude.');
end
complement=1-abs(beta_resp).^2;
if any(complement<=0)
    error('the SLR complementary spectrum is not positive.');
end

% Obtain the minimum-phase complementary response by cepstral factorisation
alpha_mag=sqrt(complement);
log_spec=fft(log(alpha_mag));
log_spec(2:fft_len/2)=2*log_spec(2:fft_len/2);
log_spec(fft_len/2+2:end)=0;
alpha_resp=exp(ifft(log_spec));
alpha_poly=fft(alpha_resp)/fft_len;
alpha=fliplr(alpha_poly(1:npts));

% Remove one Cayley-Klein rotation at a time
rf_angles=zeros(1,npts);
for slice=npts:-1:1
    if abs(alpha(slice))<=eps(norm(alpha,inf))
        error('the inverse SLR recursion encountered a singular alpha coefficient.');
    end
    ratio=beta(slice)/alpha(slice);
    cos_half=1/hypot(1,abs(ratio));
    sin_half=conj(cos_half*ratio);
    sin_size=abs(sin_half);
    if sin_size==0
        rf_angles(slice)=0;
    else
        rf_angles(slice)=2*atan2(sin_size,cos_half)*sin_half/sin_size;
    end
    if slice>1
        prev_alpha=cos_half*alpha+sin_half*beta;
        prev_beta=-conj(sin_half)*alpha+cos_half*beta;
        alpha=prev_alpha(2:slice);
        beta=prev_beta(1:(slice-1));
    end
end

% Convert per-slice rotations into Spinach control amplitudes
slice_dur=dur/npts;
rf_ctrl=rf_angles/slice_dur;
if any(~isfinite(rf_ctrl))
    error('the inverse SLR transform produced non-finite controls.');
end

% Return Cartesian, timing, and polar waveform coordinates
Cx=real(rf_ctrl);
Cy=imag(rf_ctrl);
durs=slice_dur*ones(1,npts);
amps=abs(rf_ctrl);
phis=angle(rf_ctrl);

end

% Consistency enforcement
function grumble(npts,dur,tbw,flip_angle,pass_rip,stop_rip)
if (~isnumeric(npts))||(~isreal(npts))||(~isscalar(npts))||...
   (~isfinite(npts))||(npts<2)||(mod(npts,2)~=0)
    error('npts must be a finite positive even integer.');
end
if (~isnumeric(dur))||(~isreal(dur))||(~isscalar(dur))||...
   (~isfinite(dur))||(dur<=0)
    error('dur must be a finite positive real scalar.');
end
if (~isnumeric(tbw))||(~isreal(tbw))||(~isscalar(tbw))||...
   (~isfinite(tbw))||(tbw<=0)
    error('tbw must be a finite positive real scalar.');
end
if (~isnumeric(flip_angle))||(~isreal(flip_angle))||(~isscalar(flip_angle))||...
   (~isfinite(flip_angle))||(flip_angle<=0)||(flip_angle>pi/2)
    error('flip_angle must be a finite real scalar between zero and pi/2.');
end
if (~isnumeric(pass_rip))||(~isreal(pass_rip))||(~isscalar(pass_rip))||...
   (~isfinite(pass_rip))||(pass_rip<=0)||(pass_rip>=1)
    error('pass_rip must be a finite real scalar strictly between zero and one.');
end
if (~isnumeric(stop_rip))||(~isreal(stop_rip))||(~isscalar(stop_rip))||...
   (~isfinite(stop_rip))||(stop_rip<=0)||(stop_rip>=1)
    error('stop_rip must be a finite real scalar strictly between zero and one.');
end
end

