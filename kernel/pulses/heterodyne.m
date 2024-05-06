% Signal heterodyne from wall clock time into the rotating frame. Uses
% a GPU if one is available. Syntax:
% 
%                    heterodyne(dt,exp_data,car_freq)
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
%                rotating frame
%
%
% Note: the signal must be sampled with at least four points per period
%       of the frequency being demodulated.
%
% a.acharya@soton.ac.uk
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=heterodyne.m>

function [X,Y]=heterodyne(dt,signal,freq)

% Check consistency
grumble(dt,signal,freq);

% Build time grid
time_grid=dt*((1:numel(signal))'-1);

% Move inputs to GPU
if canUseGPU()
    signal=gpuArray(signal); 
    time_grid=gpuArray(time_grid);
end

% Define a lowpass filter
F=designfilt('lowpassfir','FilterOrder',8,...
             'HalfPowerFrequency',0.25);

% Decide the output
if nargout==0

    % Plot the filter profile
    freqz(F.Coefficients,1,[],1/dt);

else

    % Mix with carrier frequency
    X=2*signal.*cos(2*pi*freq*time_grid);
    Y=2*signal.*sin(2*pi*freq*time_grid);
    clear('signal','time_grid');

    % Get coefficients
    B=F.Coefficients; 
    if canUseGPU()
        B=gpuArray(B);
    end

    % Apply the filter
    X=gather(fftfilt(B,X)); 
    Y=gather(fftfilt(B,Y));

end

end

% Consistency enforcement
function grumble(dt,signal,freq)
if (~isnumeric(dt))||(~isreal(dt))||(~isscalar(dt))||(dt<=0)
    error('dt must be a positive real number.');
end
if (~isnumeric(signal))||(~isreal(signal))||(~iscolumn(signal))
    error('signal must be a real column vector.');
end
if (~isnumeric(freq))||(~isreal(freq))||(~isscalar(freq))
    error('freq must be a real number.');
end
if 4*dt > 1/freq
    error('the specified frequency is not sampled well enough.');
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

