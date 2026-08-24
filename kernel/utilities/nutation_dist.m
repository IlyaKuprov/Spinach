% Nutation frequency distribution from a complex nutation curve. Syntax:
%
%                   [freq,distr]=nutation_dist(curve,dt)
%
% Parameters:
%
%    curve  - complex X+iY nutation curve, a row or column vector
%
%    dt     - sampling interval in seconds
%
% Outputs:
%
%    freq   - nutation frequency grid in rad/s, a column vector
%
%    distr  - non-negative nutation frequency density in inverse rad/s,
%             normalised to unit integral over freq
%
% Notes: the signal is modelled as a complex scale factor times the
%        characteristic function of a positive frequency distribution.
%        The frequency support, complex scale phase, sub-sample time
%        shift, and second-derivative Tikhonov parameter are selected
%        from the supplied trace. The returned density is therefore a
%        stable estimate rather than a raw finite-time Fourier transform.
%
%        The sign of the complex exponential is inferred from the dominant
%        Fourier half-plane, while freq is always returned as non-negative.
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=nutation_dist.m>

function [freq,distr]=nutation_dist(curve,dt)

% Check consistency
grumble(curve,dt);

% Fold the input into a column
signal=curve(:);
npts=numel(signal);
time=(0:(npts-1)).'*dt;

% Estimate the complex noise variance from the second difference
diff_two=signal(3:end)-2*signal(2:end-1)+signal(1:end-2);
noise_var=median(abs(diff_two).^2)/(6*log(2));
signal_peak=max(abs(signal).^2);
noise_var=min(noise_var,0.05*signal_peak);
noise_var=max(noise_var,eps*norm(signal,2)^2/npts);

% Apply a Hann window for support estimation
window=0.5-0.5*cos(2*pi*(0:(npts-1)).'/(npts-1));

% Build a zero-filled Fourier spectrum
nfft=2^nextpow2(8*npts);
spectrum=fftshift(fft(window.*signal,nfft));
freq_axis=(-floor(nfft/2):ceil(nfft/2)-1).'*2*pi/(nfft*dt);

% Select the Fourier half-plane containing the nutation frequencies
positive=freq_axis>0;
negative=freq_axis<0;
pos_power=sum(abs(spectrum(positive)).^2);
neg_power=sum(abs(spectrum(negative)).^2);
if pos_power>=neg_power
    phase_sign=1;
    selected=freq_axis>=0;
    selected_freq=freq_axis(selected);
    selected_power=abs(spectrum(selected)).^2;
else
    phase_sign=-1;
    selected=freq_axis<=0;
    selected_freq=-flipud(freq_axis(selected));
    selected_power=flipud(abs(spectrum(selected)).^2);
end

% Locate the noise-limited spectral support
freq_step=2*pi/(nfft*dt);
noise_bin=noise_var*sum(window.^2);
power_cut=max(10*noise_bin,1e-5*max(selected_power));
active=find(selected_power>=power_cut);
if isempty(active)
    peak_idx=find(selected_power==max(selected_power),1,'first');
    active=max(1,peak_idx-2):min(numel(selected_freq),peak_idx+2);
end
freq_lo=max(0,selected_freq(active(1))-2*freq_step);
freq_hi=min(pi/dt,selected_freq(active(end))+2*freq_step);
if freq_hi<=freq_lo
    freq_lo=0;
    freq_hi=min(pi/dt,selected_freq(active(end))+4*freq_step);
end
if freq_hi<=freq_lo
    error('the curve does not contain a resolvable frequency range.');
end

% Choose a frequency grid at approximately twice the Fourier resolution
resolution=ceil((freq_hi-freq_lo)*npts*dt/(2*pi));
ngrid=min(400,max(160,2*resolution));
freq=linspace(freq_lo,freq_hi,ngrid).';

% Build the second-difference regularisation matrix
D=spdiags(ones(ngrid-2,1)*[1 -2 1],0:2,ngrid-2,ngrid);

% Scale the regularisation search to the data and grid
kernel=exp(phase_sign*1i*time*freq.');
real_kernel=[real(kernel);imag(kernel)];
scale=trace(real_kernel'*real_kernel)/max(trace(D'*D),eps);
lambdas=scale*logspace(-8,4,13);

% Choose the strongest regularisation within 10% of the best misfit
trial_err=zeros(size(lambdas));
for n=1:numel(lambdas)
    [~,trial_err(n)]=fit_nutation(signal,kernel,real_kernel,D,lambdas(n));
end
lambda=lambdas(find(trial_err<=1.1*min(trial_err),1,'last'));

% Search for a small error in the time origin
shift_grid=linspace(-2*dt,2*dt,9);
best_err=inf;
best_shift=0;
for n=1:numel(shift_grid)
    shifted_time=time+shift_grid(n);
    kernel=exp(phase_sign*1i*shifted_time*freq.');
    real_kernel=[real(kernel);imag(kernel)];
    [~,fit_err]=fit_nutation(signal,kernel,real_kernel,D,lambda);
    if fit_err<best_err
        best_err=fit_err;
        best_shift=shift_grid(n);
    end
end

% Refine the best time origin on the local interval
shift_step=shift_grid(2)-shift_grid(1);
shift_grid=linspace(best_shift-shift_step,best_shift+shift_step,5);
for n=1:numel(shift_grid)
    shifted_time=time+shift_grid(n);
    kernel=exp(phase_sign*1i*shifted_time*freq.');
    real_kernel=[real(kernel);imag(kernel)];
    [~,fit_err]=fit_nutation(signal,kernel,real_kernel,D,lambda);
    if fit_err<best_err
        best_err=fit_err;
        best_shift=shift_grid(n);
    end
end

% Re-select the regularisation at the corrected time origin
shifted_time=time+best_shift;
kernel=exp(phase_sign*1i*shifted_time*freq.');
real_kernel=[real(kernel);imag(kernel)];
for n=1:numel(lambdas)
    [~,trial_err(n)]=fit_nutation(signal,kernel,real_kernel,D,lambdas(n));
end
lambda=lambdas(find(trial_err<=1.1*min(trial_err),1,'last'));

% Final fit at the corrected origin and regularisation
best_weights=fit_nutation(signal,kernel,real_kernel,D,lambda);

% Convert fitted probability masses into a unit-integral density
distr=max(best_weights,0);
distr=distr/trapz(freq,distr);

end

% Fit a non-negative regularised distribution with an unknown complex scale
function [weights,fit_err]=fit_nutation(signal,kernel,real_kernel,D,lambda)

% Stack the data fit and the regularisation penalty
lsq_mat=[real_kernel;sqrt(lambda)*full(D)];
zero_pen=zeros(size(D,1),1);
phase=exp(1i*angle(signal(1)));

% Alternate phase estimation and non-negative Tikhonov fitting
for n=1:8
    rotated=conj(phase)*signal;
    weights=lsqnonneg(lsq_mat,[real(rotated);imag(rotated);zero_pen]);
    predicted=kernel*weights;
    overlap=predicted'*signal;
    if abs(overlap)>0
        phase=overlap/abs(overlap);
    else
        phase=1;
    end
end

% Evaluate the fitted residual
fit_err=norm(phase*predicted-signal,2)^2;

end

% Consistency enforcement
function grumble(curve,dt)
if (~isnumeric(curve))||(~isvector(curve))||isempty(curve)||isreal(curve)
    error('curve must be a non-empty complex numeric vector.');
end
if numel(curve)<8
    error('curve must contain at least eight points.');
end
if any(~isfinite(curve))
    error('curve must not contain Inf or NaN.');
end
if (~isnumeric(dt))||(~isreal(dt))||(~isscalar(dt))||(~isfinite(dt))||(dt<=0)
    error('dt must be a positive real scalar.');
end
end


