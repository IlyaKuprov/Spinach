% Frequency-swept inversion pulse with offset-independent adiabati-
% city (Tannus and Garwood, JMR A 120, 133 (1996)). The user supp-
% lies an amplitude modulation function; the frequency sweep that
% makes the adiabaticity factor the same for every offset inside
% the sweep bandwidth is obtained by integrating the square of that
% amplitude function. Syntax:
%
%    [Cx,Cy,durs,ints,amps,phis,frqs]=...
%                      oia_pulse(npts,dur,bwidth,am_fun)
%
% Parameters:
%
%        npts    - number of discretisation points in
%                  the waveform, a floating-point sca-
%                  lar integer not smaller than 2
%
%        dur     - pulse duration, seconds, a positive
%                  floating-point scalar
%
%      bwidth    - sweep bandwidth around zero frequ-
%                  ency, Hz, a positive floating-point
%                  scalar
%
%      am_fun    - amplitude modulation function handle,
%                  F1(tau) in Table 1 of the paper, that
%                  accepts a row vector of normalised ti-
%                  mes in the [-1,1] interval and returns
%                  a row vector of non-negative floating-
%                  point amplitudes that are positive at
%                  all interior points, the scale of which
%                  does not matter
%                  because the envelope is normalised to
%                  unit peak amplitude internally
%
% Outputs:
%
%          Cx    - real part of the waveform, calibrated to
%                  the same adiabaticity factor as the
%                  inversion pulse in chirp_pulse.m, rad/s,
%                  a row vector with npts elements
%
%          Cy    - imag part of the waveform, calibrated to
%                  the same adiabaticity factor as the
%                  inversion pulse in chirp_pulse.m, rad/s,
%                  a row vector with npts elements
%
%        durs    - slice durations for piecewise-constant
%                  approximation, seconds, a row vector
%                  with npts elements
%
%        ints    - interval durations for piecewise-linear
%                  approximation, seconds, a row vector
%                  with npts-1 elements
%
%        amps    - waveform amplitudes, rad/s, a row
%                  vector with npts elements
%
%        phis    - waveform phases, radians, zero at the
%                  centre of the pulse, a row vector with
%                  npts elements
%
%        frqs    - waveform frequencies, Hz, a row vector
%                  with npts elements
%
% Note: the amplitude functions in Table 1 of the paper, written
%       with the 1% edge truncation used there, are
%
%           Lorentz    @(tau)1./(1+99*tau.^2)
%           HS         @(tau)sech(asech(0.01)*tau)
%           Gauss      @(tau)exp(-log(100)*tau.^2)
%           Hanning    @(tau)(1+cos(pi*tau))/2
%           HSn        @(tau)sech(asech(0.01)*tau.^n)
%           Sin^n      @(tau)1-abs(sin(pi*tau/2)).^n
%
%       and the hand drawn pulse in Figure 3 of the paper is any
%       smooth function that is positive inside the pulse and
%       may vanish only at its two ends.
%
% Note: the frequency sweep runs from -bwidth/2 to +bwidth/2; for
%       a chirp pulse, which has a constant amplitude function, the
%       output of this function coincides with chirp_pulse.m.
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=oia_pulse.m>

function [Cx,Cy,durs,ints,amps,phis,frqs]=oia_pulse(npts,dur,bwidth,am_fun)

% Check consistency
grumble(npts,dur,bwidth,am_fun);

% Normalised time grid
time_grid=linspace(-0.5,0.5,npts);
ints=dur*diff(time_grid);
durs=(dur/npts)*ones(1,npts);

% Amplitude function on the [-1,1] interval, unit peak
am_vals=am_fun(2*time_grid);
if (~isfloat(am_vals))||(~isreal(am_vals))||(~isrow(am_vals))||...
   (numel(am_vals)~=npts)||(~all(isfinite(am_vals)))||any(am_vals<0)||...
   all(am_vals==0)||any(am_vals(2:end-1)==0)
    error('am_fun must return a row of non-negative finite floating-point numbers, positive inside the pulse.');
end
am_vals=am_vals/max(am_vals);

% Sweep rate proportional to the amplitude squared, Eq. 6 in the paper
sweep_cdf=cumtrapz(time_grid,am_vals.^2);
sweep_cdf=sweep_cdf/sweep_cdf(end);

% Frequencies and phases, zero phase at the centre of the pulse
frqs=bwidth*(sweep_cdf-0.5);
phis=2*pi*dur*cumtrapz(time_grid,frqs);
phis=phis-interp1(time_grid,phis,0,'spline');

% Check sampling adequacy
phi_jumps=2*pi*max(abs(frqs(1:end-1)),abs(frqs(2:end))).*ints;
if any(abs(phi_jumps)>pi,'all')&&(nargout<7)
    error('insufficient number of points to sample the pulse.');
end

% Amplitude with the same adiabaticity factor as chirp_pulse.m
amps=2*pi*sqrt(bwidth/(dur*trapz(time_grid,am_vals.^2)))*am_vals;

% Convert into Cartesian coordinates
[Cx,Cy]=polar2cartesian(amps,phis);

end

% Consistency enforcement
function grumble(npts,dur,bwidth,am_fun)
if (~isfloat(npts))||(~isreal(npts))||(numel(npts)~=1)||...
   (~isfinite(npts))||(npts<2)||(mod(npts,1)~=0)
    error('npts must be a finite real floating-point integer greater than 1.');
end
if (~isfloat(dur))||(~isreal(dur))||...
   (numel(dur)~=1)||(~isfinite(dur))||(dur<=0)
    error('dur must be a finite positive real floating-point number.');
end
if (~isfloat(bwidth))||(~isreal(bwidth))||...
   (numel(bwidth)~=1)||(~isfinite(bwidth))||(bwidth<=0)
    error('bwidth must be a finite positive real floating-point number.');
end
if ~isa(am_fun,'function_handle')
    error('am_fun must be a function handle.');
end
end

% One of the authors (A.T.) asked his 11-year-old son to hand
% draw an arbitrary F1(t) function. He was instructed to start
% and end at amplitude zero and to keep the shape relatively
% smooth avoiding steep changes.
%
% Alberto Tannus and Michael Garwood,
% J. Magn. Reson. A 120, 133 (1996)


