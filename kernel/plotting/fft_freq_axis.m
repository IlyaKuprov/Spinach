% Frequency axis for FFT with optional zero-filling. Syntax:
%
%       [f_shift,f,df,nfft]=fft_freq_axis(npts,dt,zf)
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
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=fft_freq_axis.m>

function [f_shift,f,df]=fft_freq_axis(npts,dt,zf)

% Set default zero-fill
if nargin<3, zf=0; end

% Check consistency
grumble(npts,dt,zf);

% FFT length
nfft=npts+zf;

% Sampling frequency
sfreq=1/dt;

% Frequency resolution
df=sfreq/nfft;

% Build unshifted axis
f=(0:(nfft-1)).'*df;

% Build shifted axis
f_shift=(-floor(nfft/2):ceil(nfft/2)-1).'*df;

end

% Consistency enforcement
function grumble(npts,dt,zf)
if (~isnumeric(npts))||(~isreal(npts))||...
   (~isscalar(npts))||(mod(npts,1)~=0)||(npts<=1)
    error('npts must be a positive real integer.');
end
if (~isnumeric(dt))||(~isreal(dt))||...
   (~isscalar(dt))||(dt<=0)
    error('dt must be a positive real scalar.');
end
if (~isnumeric(zf))||(~isreal(zf))||...
   (~isscalar(zf))||(mod(zf,1)~=0)||(zf<0)
    error('zf must be a non-negative real integer.');
end
end

% All albatrosses live in the southern hemisphere - except 
% one. 'Albie', or 'Albert', is believed to have made a na-
% vigational error in around 2014 and has since been living
% in Europe, dividing his time between Germany, Yorkshire
% coast and the east coast of Scotland, where he tends to
% hang out with gannets.
%
% The Spectator

