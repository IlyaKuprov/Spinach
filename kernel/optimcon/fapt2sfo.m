% Converts a freq-ampl-phase-time specification of a pulse sequ-
% uence into the corresponding single frequency origin waveform
% that is compatible with GRAPE optimisations. Syntax:
%
%         [wave,dt,time_grid]=fapt2sfo(fapt,time_grid)
%
% Parameters:
%
%    fapt  - a cell array of 5-element row vectors with of 
%            the following structure: [frequency (Hz), amp-
%            litude (rad/s), phase at t=0 (radians), start
%            time (seconds), end time (seconds)]
%
%    time_grid - optional vector of time grid ticks; when
%                not provided, the grid is made at twice 
%                the Nyquist-Shannon minimum sampling ra-
%                te of the highest frequency present
%
% Output:
%
%    wave  - pulse sequence as a single waveform, a matrix
%            with two rows, corresponding to X and Y compo-
%            nents
%
%    dt    - step duration of the time grid, seconds
%
%    time_grid - row vector of time grid ticks
%
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=fapt2sfo.m>

function [wave,dt,time_grid]=fapt2sfo(fapt,time_grid)

% Check consistency
grumble(fapt);

% Make time grid if not provided
if ~exist('time_grid','var')

    % Find the highest frequency in the sequence
    max_freq=max(abs(cellfun(@(x)x(1),fapt)));

    % Find the end time of the sequence
    end_time=max(cellfun(@(x)x(5),fapt));

    % Get the time grid with dt < half Nyquist
    dt_hN=1/(4*max_freq); npts=ceil(end_time/dt_hN)+1;
    time_grid=linspace(0,end_time,npts); dt=time_grid(2);

else

    % Accept what is given
    npts=numel(time_grid); dt=[];

end

% Preallocate waveform
wave=zeros(2,npts);

% Loop over events
for n=1:numel(fapt)
    
    % Pull the parameters out
    freq=fapt{n}(1); ampl=fapt{n}(2); phi=fapt{n}(3);
    start_time=fapt{n}(4); end_time=fapt{n}(5);
    
    % Get the time interval mask
    time_mask=(time_grid>=start_time)&(time_grid<=end_time);

    % Add event to the waveform, counterclockwise phase
    wave(1,:)=wave(1,:)+ampl*cos(2*pi*freq*time_grid+phi).*time_mask;
    wave(2,:)=wave(2,:)+ampl*sin(2*pi*freq*time_grid+phi).*time_mask;

end

end

% Consistency enforcement
function grumble(fapt)
if ~iscell(fapt)
    error('fapt must be a cell array of five-element vectors.');
end
for n=1:numel(fapt)
    if (~isnumeric(fapt{n}))||(~isreal(fapt{n}))||(numel(fapt{n})~=5)
        error('fapt must be a cell array of real five-element vectors.');
    end
    if fapt{n}(2)<0
        error('negative amplitude found in fapt.');
    end
    if fapt{n}(4)>=fapt{n}(5)
        error('end time precedes start time in fapt.');
    end
end
end

% With four parameters I can fit an elephant, and with 
% five I can make him wiggle his trunk.
%
% John von Neumann

