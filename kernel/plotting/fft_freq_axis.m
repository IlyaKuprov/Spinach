% Frequency axis for FFT with optional zero-filling. Syntax:
%
%          [f_shift,f,df,nfft]=fft_freq_axis(npts,dt,zf)
%
% Parameters:
%
%    npts    - number of acquired time-domain points
%
%    dt      - time step between points
%
%    zf      - zero-fill length added to time domain
%
% Outputs:
%
%    f_shift - frequency axis for fftshift(fft(...))
%
%    f       - frequency axis for fft(...)
%
%    df      - frequency resolution
%
%    nfft    - FFT length used

function [f_shift,f,df,nfft]=fft_freq_axis(npts,dt,zf)

% Check consistency
grumble(npts,dt,zf);

% Set default zero-fill
if nargin<3, zf=0; end

% Round point counts
npts=round(npts);
zf=round(zf);

% Compute FFT length
nfft=npts+zf;

% Compute sampling frequency
sfreq=1/dt;

% Compute frequency resolution
df=sfreq/nfft;

% Build unshifted axis
f=(0:(nfft-1)).'*df;

% Build shifted axis
f_shift=(-floor(nfft/2):ceil(nfft/2)-1).'*df;

end

% Consistency enforcement
function grumble(npts,dt,zf)
if (~isnumeric(npts))||(~isreal(npts))||...
   (~isscalar(npts))||(npts<=0)
    error('npts must be a positive real scalar.');
end
if (~isnumeric(dt))||(~isreal(dt))||...
   (~isscalar(dt))||(dt<=0)
    error('dt must be a positive real scalar.');
end
if nargin<3
    return;
end
if (~isnumeric(zf))||(~isreal(zf))||...
   (~isscalar(zf))||(zf<0)
    error('zf must be a non-negative real scalar.');
end
end
