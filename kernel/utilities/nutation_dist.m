% Nutation frequency distribution from a nutation curve measured with
% the same coil used for excitation and detection. Syntax:
%
%                [freq,distr]=nutation_dist(curve,dt,lambda)
%
% Parameters:
%
%    curve  - nutation curve, a row or column vector; either
%             the complex X+iY output of a quadrature receiver,
%             or a real phase-corrected trace
%
%    dt     - sampling interval in seconds
%
%    lambda - second-derivative Tikhonov regularisation para-
%             meter, a non-negative real scalar; the curve is
%             normalised to unit maximum modulus and the fit-
%             ting kernel is dimensionless, so lambda is a di-
%             mensionless number of the order of the ratio of
%             the squared Frobenius norms of the kernel and of
%             the second difference matrix; zero switches the
%             regularisation off
%
% Outputs:
%
%    freq   - nutation frequency grid in rad/s, a column vector
%
%    distr  - non-negative nutation frequency density in inverse
%             rad/s, normalised to unit integral over freq
%
% Notes: when one coil both excites and detects, the reciprocity
%        principle makes the detected amplitude of every isochromat
%        proportional to its own nutation frequency, and only the
%        sine component of the nutation is observable. The curve is
%        therefore modelled as an unknown complex receiver scale
%        times
%
%           s(t)=integral(distr(freq)*freq*sin(freq*(t+t0)),d freq)
%
%        and the reception weight is divided out, so that the
%        returned density is the true nutation frequency distribu-
%        tion. A real receiver scale is a valid special case, and
%        a phase-corrected real trace is therefore accepted. The
%        frequency support, receiver phase, and sub-sample time
%        shift are selected from the supplied trace; the returned
%        density is therefore a stable estimate rather than a raw
%        finite-time transform. The reconstruction runs on twice
%        the noise-limited support width, centred on the same band
%        and clipped to the Nyquist interval, so that the margins
%        of the distribution are visible rather than truncated.
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=nutation_dist.m>

function [freq,distr]=nutation_dist(curve,dt,lambda)

% Check consistency
grumble(curve,dt,lambda);

% Fold the input into a column
signal=curve(:);
npts=numel(signal);
time=(0:(npts-1)).'*dt;

% Normalise the receiver gain away
signal=signal/norm(signal,inf);

% Estimate the complex noise variance from the second difference
diff_two=signal(3:end)-2*signal(2:end-1)+signal(1:end-2);
noise_var=median(abs(diff_two).^2)/(6*log(2));
signal_peak=max(abs(signal).^2);
noise_var=min(noise_var,0.05*signal_peak);
noise_var=max(noise_var,eps*norm(signal,2)^2/npts);

% Rotate the receiver phase line onto the real axis
line_phase=exp(1i*angle(sum(signal.^2))/2);
rotated=real(conj(line_phase)*signal);

% Apply a fade-out window for support estimation
window=0.5+0.5*cos(pi*(0:(npts-1)).'/(npts-1));

% Sine-transform power spectrum on a zero-filled grid
nfft=2^nextpow2(8*npts);
spectrum=fftshift(fft(window.*rotated,nfft));
freq_axis=(-floor(nfft/2):ceil(nfft/2)-1).'*2*pi/(nfft*dt);
selected=freq_axis>=0;
selected_freq=freq_axis(selected);
selected_power=imag(spectrum(selected)).^2;

% Locate the noise-limited spectral support
freq_step=2*pi/(nfft*dt);
noise_bin=noise_var*sum(window.^2)/4;
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

% Double the support about its centre to expose the distribution margins
freq_mid=(freq_lo+freq_hi)/2;
freq_span=freq_hi-freq_lo;
freq_lo=max(0,freq_mid-freq_span);
freq_hi=min(pi/dt,freq_mid+freq_span);

% Choose a frequency grid at approximately twice the Fourier resolution
resolution=ceil((freq_hi-freq_lo)*npts*dt/(2*pi));
ngrid=min(400,max(160,2*resolution));
freq=linspace(freq_lo,freq_hi,ngrid).';

% Second-difference regularisation matrix with clamped boundaries
D=spdiags(ones(ngrid,1)*[1 -2 1],-1:1,ngrid,ngrid);

% Dimensionless reciprocity reception weight
rec_weight=freq.'/freq_hi;

% Silence the non-negative least squares solver
options=optimset('Display','off');

% Search for a small error in the time origin
shift_grid=linspace(-2*dt,2*dt,9);
best_err=inf;
best_shift=0;
for n=1:numel(shift_grid)
    kernel=sin((time+shift_grid(n))*freq.').*rec_weight;
    [~,fit_err]=fit_nutation(signal,kernel,D,lambda,options);
    if fit_err<best_err
        best_err=fit_err;
        best_shift=shift_grid(n);
    end
end

% Refine the best time origin on the local interval
shift_step=shift_grid(2)-shift_grid(1);
shift_grid=linspace(best_shift-shift_step,best_shift+shift_step,5);
for n=1:numel(shift_grid)
    kernel=sin((time+shift_grid(n))*freq.').*rec_weight;
    [~,fit_err]=fit_nutation(signal,kernel,D,lambda,options);
    if fit_err<best_err
        best_err=fit_err;
        best_shift=shift_grid(n);
    end
end

% Final fit at the corrected time origin
kernel=sin((time+best_shift)*freq.').*rec_weight;
best_weights=fit_nutation(signal,kernel,D,lambda,options);

% Convert fitted probability masses into a unit-integral density
distr=max(best_weights,0);
distr=distr/trapz(freq,distr);

end

% Fit a non-negative regularised distribution with an unknown receiver phase
function [weights,fit_err]=fit_nutation(signal,kernel,D,lambda,options)

% Stack the data fit and the regularisation penalty
lsq_mat=[kernel;sqrt(lambda)*full(D)];
zero_pen=zeros(size(D,1),1);

% Estimate the receiver phase line with a pi ambiguity
line_phase=exp(1i*angle(sum(signal.^2))/2);

% Try both signs of the receiver phase estimate
fit_err=inf;
for candidate=[line_phase -line_phase]

    % Alternate phase estimation and non-negative Tikhonov fitting
    phase=candidate;
    for n=1:4
        rotated=real(conj(phase)*signal);
        trial_weights=lsqnonneg(lsq_mat,[rotated;zero_pen],options);
        predicted=kernel*trial_weights;
        overlap=predicted'*signal;
        if abs(overlap)>0
            phase=overlap/abs(overlap);
        else
            phase=candidate;
        end
    end

    % Keep the better of the two phase branches
    trial_err=norm(phase*predicted-signal,2)^2;
    if trial_err<fit_err
        fit_err=trial_err; weights=trial_weights;
    end

end

end

% Consistency enforcement
function grumble(curve,dt,lambda)
if (~isnumeric(curve))||(~isvector(curve))||isempty(curve)
    error('curve must be a non-empty numeric vector.');
end
if numel(curve)<8
    error('curve must contain at least eight points.');
end
if any(~isfinite(curve))
    error('curve must not contain Inf or NaN.');
end
if ~any(curve)
    error('curve must not be identically zero.');
end
if (~isnumeric(dt))||(~isreal(dt))||(~isscalar(dt))||(~isfinite(dt))||(dt<=0)
    error('dt must be a positive real scalar.');
end
if (~isnumeric(lambda))||(~isreal(lambda))||(~isscalar(lambda))||...
   (~isfinite(lambda))||(lambda<0)
    error('lambda must be a non-negative real scalar.');
end
end


