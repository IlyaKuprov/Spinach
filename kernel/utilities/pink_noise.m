% Trajectories of a real random process with a 1/f power spectral
% density S(omega)=amp^2/|omega/2pi| between an infrared and an ul-
% traviolet cutoff, synthesised in the frequency domain: every Fou-
% rier component within the band receives a complex Gaussian ampli-
% tude with the variance set by the spectral density, and the inver-
% se FFT returns the time domain trajectory (Timmer and Koenig, As-
% tron. Astrophys. 300, 707, 1995). Syntax:
%
%           trajs=pink_noise(amp,dt,npts,ntraj,f_ir,f_uv)
%
% Parameters:
%
%    amp    - noise amplitude, the two-sided spectral density
%             is amp^2/|f| in units of [amp]^2/Hz
%
%    dt     - sampling interval, seconds
%
%    npts   - number of points in each trajectory
%
%    ntraj  - number of trajectories
%
%    f_ir   - infrared cutoff frequency, Hz, not below the
%             frequency resolution 1/(npts*dt)
%
%    f_uv   - ultraviolet cutoff frequency, Hz, not above the
%             Nyquist frequency 1/(2*dt)
%
% Outputs:
%
%    trajs  - trajectories, [ntraj x npts] real array, the
%             variance of each point is 2*amp^2*log(f_uv/f_ir)
%
% Note: the trajectories are periodic with the period npts*dt;
%       use a duration well above 1/f_ir to avoid artefacts.
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=pink_noise.m>

function trajs=pink_noise(amp,dt,npts,ntraj,f_ir,f_uv)

% Check consistency
grumble(amp,dt,npts,ntraj,f_ir,f_uv);

% Positive frequency bins of the grid
nbins=floor(npts/2); freqs=(1:nbins)/(npts*dt);

% Standard deviation of the quadrature amplitudes in the band
sigma=zeros(1,nbins); band=(freqs>=f_ir)&(freqs<=f_uv);
sigma(band)=sqrt(2*amp^2./(freqs(band)*npts*dt));

% Complex Gaussian spectrum in the MATLAB ifft normalisation
spec=zeros(ntraj,npts);
spec(:,2:(nbins+1))=(npts/2)*(randn(ntraj,nbins)+1i*randn(ntraj,nbins)).*sigma;

% Real trajectories from the positive frequency half
trajs=2*real(ifft(spec,[],2));

end

% Consistency enforcement
function grumble(amp,dt,npts,ntraj,f_ir,f_uv)
if (~isnumeric(amp))||(~isreal(amp))||(~isscalar(amp))||(amp<=0)
    error('amp must be a positive real scalar.');
end
if (~isnumeric(dt))||(~isreal(dt))||(~isscalar(dt))||(dt<=0)
    error('dt must be a positive real scalar.');
end
if (~isnumeric(npts))||(~isreal(npts))||(~isscalar(npts))||(npts<4)||(mod(npts,1)~=0)
    error('npts must be an integer greater than 3.');
end
if (~isnumeric(ntraj))||(~isreal(ntraj))||(~isscalar(ntraj))||(ntraj<1)||(mod(ntraj,1)~=0)
    error('ntraj must be a positive integer.');
end
if (~isnumeric(f_ir))||(~isreal(f_ir))||(~isscalar(f_ir))||(f_ir<1/(npts*dt))
    error('f_ir must be a real scalar not below 1/(npts*dt).');
end
if (~isnumeric(f_uv))||(~isreal(f_uv))||(~isscalar(f_uv))||(f_uv<=f_ir)||(f_uv>1/(2*dt))
    error('f_uv must be a real scalar between f_ir and 1/(2*dt).');
end
end

% The noise is the signal.
%
% (a proverb of every spectroscopist who has
%  ever tried to measure a relaxation time)

