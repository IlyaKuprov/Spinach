% Hyperbolic secant pulse in Cartesian and amplitude-phase
% representation. Syntax:
%
%    [Cx,Cy,time_grid,amps,phis]=...
%     sech_pulse(peak_ampl,freq_mod,phase_mod,dur,npts)
%
% Parameters:
%
%     peak_amp  - peak amplitude, rad/s
%
%     freq_mod  - frequency modulation parameter, rad/s
%
%     phase_mod - phase modulation parameter, dimless
%
%     dur       - pulse duration, seconds
%
%     npts      - number of digitisation points
%
% Outputs:
%
%     Cx        - a vector of coefficients in front of Sx
%                 spin operator at each time slice, rad/s
%
%     Cy        - a vector of coefficients in front of Sy
%                 spin operator at each time slice, rad/s
%
%     time_grid - a vector of time grid points, seconds
%
%     amps      - a vector of pulse amplitudes at each ti-
%                 me slice, rad/s
%
%     phis      - a vector of pulse phases at each time 
%                 slice (phi=0 at t=0 in the centre), rad
%
% Example: 
%
%     [Cx,Cy,time_grid]=sech_pulse(1,672,5,10.24e-3,1000);
%     plot(time_grid,[Cx; Cy]); kgrid; xlim tight;
%     kxlabel('time, seconds'); kylabel('amplitude, rad/s');
%
% ilya.kuprov@weizmann.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=sech_pulse.m>

function [Cx,Cy,time_grid,amps,phis]=...
          sech_pulse(peak_amp,freq_mod,phase_mod,dur,npts)

% Check consistency
grumble(peak_amp,freq_mod,phase_mod,dur,npts);

% Get the time grid, t=0 is centre
time_grid=linspace(-dur/2,dur/2,npts);

% Get the amplitudes
amps=peak_amp*sech(freq_mod*time_grid);

% Get the phases (phi=0 at t=0)
phis=phase_mod*log(cosh(freq_mod*time_grid));

% Convert into Cartesians
[Cx,Cy]=polar2cartesian(amps,phis);

end

% Consistency enforcement
function grumble(peak_amp,freq_mod,phase_mod,dur,npts)
if (~isnumeric(peak_amp))||(~isreal(peak_amp))||(~isscalar(peak_amp))
    error('peak_amp must be a real scalar.');
end
if (~isnumeric(freq_mod))||(~isreal(freq_mod))||(~isscalar(freq_mod))
    error('freq_mod must be a real scalar.');
end
if (~isnumeric(phase_mod))||(~isreal(phase_mod))||(~isscalar(phase_mod))
    error('phase_mod must be a real scalar.');
end
if (~isnumeric(dur))||(~isreal(dur))||(~isscalar(dur))||(dur<=0)
    error('dur must be a positive real scalar.');
end
if (~isnumeric(npts))||(~isreal(npts))||(~isscalar(npts))||...
   (npts<1)||(mod(npts,1)~=0)
    error('npts must be a positive integer.');
end

end

% We first sent it to a Nature journal, and within 24 hours 
% they rejected it as an incremental contribution. I started
% learning English only at university, so I had to look up
% the meaning of the word 'incremental'.
%
% Katalin Kaliko, about her paper that
% got her the Nobel Prize in 2023

