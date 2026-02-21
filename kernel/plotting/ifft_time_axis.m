% Time axis for IFFT with optional zero-filling. Syntax:
%
%       [t_shift,t,dt,nifft]=ifft_time_axis(npts,df,zf)
%
% Parameters:
%
%    npts    - number of frequency-domain points
%
%    df      - frequency interval between points, Hz
%
%    zf      - zero-fill length added to either 
%              side of the frequency domain
%
% Outputs:
%
%    t_shift - time axis for fftshift(ifft(...))
%
%    t       - time axis for ifft(...)
%
%    dt      - time step between points
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=ifft_time_axis.m>

function [t_shift,t,dt]=ifft_time_axis(npts,df,zf)

% Set default zero-fill
if nargin<3, zf=0; end

% Check consistency
grumble(npts,df,zf);

% IFFT length
nifft=npts+2*zf;

% Time resolution
dt=1/(df*nifft);

% Build unshifted axis
t=(0:(nifft-1)).'*dt;

% Build shifted axis
t_shift=(-floor(nifft/2):ceil(nifft/2)-1).'*dt;

end

% Consistency enforcement
function grumble(npts,df,zf)
if (~isnumeric(npts))||(~isreal(npts))||...
   (~isscalar(npts))||(mod(npts,1)~=0)||(npts<=1)
    error('npts must be a positive real integer.');
end
if (~isnumeric(df))||(~isreal(df))||...
   (~isscalar(df))||(df<=0)
    error('df must be a positive real scalar.');
end
if (~isnumeric(zf))||(~isreal(zf))||...
   (~isscalar(zf))||(mod(zf,1)~=0)||(zf<0)
    error('zf must be a non-negative real integer.');
end
end

% Deeply amused by those telling me I've lost
% their admiration due to the disrespect I
% show [...]. I shall file your lost admira-
% tion carefully in the box where I keep my
% missing fucks.
%
% J.K. Rowling

