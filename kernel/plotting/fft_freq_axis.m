%FFT_FREQ_AXIS  Frequency axis for FFT with optional zero-filling.
%
% Inputs
%   Nt : number of acquired time-domain points
%   dt : time spacing between points (seconds, or whatever your time unit is)
%   zf : extent of zero-filling = number of zeros appended in time domain
%        (so FFT length is Nfft = Nt + zf)
%
% Outputs
%   fShift : frequency axis (same length as FFT) for fftshift(fft(...))
%   f      : frequency axis for raw fft(...) (unshifted, starts at 0)
%   df     : frequency resolution
%   Nfft   : FFT length used

function [fShift, f, df, Nfft] = fft_freq_axis(Nt, dt, zf)

if nargin < 3, zf = 0; end
Nt  = round(Nt);
zf  = round(zf);
if Nt <= 0, error('Nt must be positive.'); end
if dt <= 0, error('dt must be positive.'); end
if zf < 0, error('zf must be >= 0.'); end

Nfft = Nt + zf;

Fs = 1/dt;        % sampling frequency (units: 1/time)
df = Fs/Nfft;     % frequency bin spacing

% Unshifted frequency bins (matches fft output ordering)
f = (0:Nfft-1).' * df;

% Shifted frequency bins (matches fftshift(fft(...)) ordering)
fShift = (-floor(Nfft/2) : ceil(Nfft/2)-1).' * df;

end
