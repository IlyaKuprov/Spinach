% Signal heterodyne from wall clock time into the rotating frame
% using analytic signal demodulation: the negative frequency half
% of the spectrum (the counter-rotating component), as well as the
% direction-ambiguous DC and Nyquist bins, are dropped; the posi-
% tive frequency half is doubled and frequency-shifted. Syntax:
%
%               [X,Y]=heterodyne(dt,signal,freq)
%
% Parameters:
%
%   dt         - time step in the input data, seconds
%
%   signal     - wall clock time signal, a column vector
%
%   freq       - frequency to be demodulated, Hz
%
% Outputs:
%
%   X, Y       - in-phase and out-of-phase parts of the
%                rotating frame signal, column vectors
%
% Notes: the signal must be sampled with more than two points per
%        period of the frequency being demodulated; the transform
%        is zero-phase, so sample k of the outputs refers to the
%        same wall clock time as sample k of the input; dropping
%        the DC bin subtracts the signal mean exactly; the record
%        is treated as periodic by the FFT, and should therefore
%        begin and end in dead time.
%
% a.acharya@soton.ac.uk
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=heterodyne.m>

function [X,Y]=heterodyne(dt,signal,freq)

% Check consistency
grumble(dt,signal,freq);

% Build time grid
time_grid=dt*((1:numel(signal))'-1);

% One-sided spectral mask, DC and Nyquist bins dropped
mask=zeros(numel(signal),1);
mask(2:ceil(numel(signal)/2))=2;

% Demodulate the analytic signal into the rotating frame
signal=ifft(mask.*fft(signal)).*exp(-2i*pi*freq*time_grid);

% In-phase and out-of-phase components
X=real(signal); Y=-imag(signal);

end

% Consistency enforcement
function grumble(dt,signal,freq)
if (~isnumeric(dt))||(~isreal(dt))||(~isscalar(dt))||...
   (~isfinite(dt))||(dt<=0)
    error('dt must be a finite positive real number.');
end
if (~isnumeric(signal))||(~isreal(signal))||(~iscolumn(signal))
    error('signal must be a real column vector.');
end
if (~isnumeric(freq))||(~isreal(freq))||(~isscalar(freq))||(~isfinite(freq))
    error('freq must be a finite real number.');
end
if (freq~=0)&&((2*dt)>=(1/abs(freq)))
    error('the carrier must be sampled with more than two points per period.');
end
end

% How did we spend the entire day monitoring what we thought was an NMR
% signal before noticing we had no field? Probably, in the slightly pa-
% nicked atmosphere a small radio station spike was mistaken for an NMR
% signal. Under normal operation with a cold superconducting magnet the
% magnet coil acts as an excellent radio frequency shield. When the mag-
% net quenches and warms up, the shielding effect is removed, so radio
% station spikes can appear. That's what seems to have misled us. Moral:
% check that the "NMR signal" disappears when no pulse is applied.
%
% Malcolm Levitt's email to
% his group, July 2023.

