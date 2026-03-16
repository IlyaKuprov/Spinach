% Returns a saw-tooth waveform. Syntax:
%
%      waveform=sawtooth(amplitude,frequency,time_grid)
%
% Arguments:
%
%      amplitude  -  amplitude at the tooth top
%
%      frequency  -  waveform frequency in teeth per second
%
%      time_grid  -  grid of time points, seconds
%
% Outputs:
%
%      waveform   -  waveform array of the same shape
%                    as the time_grid input
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=sawtooth.m>

function waveform=sawtooth(amplitude,frequency,time_grid)

% Check consistency
grumble(amplitude,frequency,time_grid);

% Compute the waveform
waveform=amplitude*(2*frequency*mod(time_grid,1/frequency)-1);

end

% Consistency enforcement
function grumble(amplitude,frequency,time_grid)
if (numel(amplitude)~=1)||(~isnumeric(amplitude))||...
   (~isreal(amplitude))||(~isfinite(amplitude))
    error('amplitude parameter must be a finite real scalar.');
end
if (numel(frequency)~=1)||(~isnumeric(frequency))||...
   (~isreal(frequency))||(~isfinite(frequency))||(frequency<=0)
    error('frequency parameter must be a positive finite real scalar.');
end
if (~isnumeric(time_grid))||(~isreal(time_grid))||...
   any(~isfinite(time_grid),'all')
    error('time_grid parameter must be a finite real array.');
end
end

% They're all aristocrats, that's true... because they know that there's
% no such thing as a lousy job - only lousy men who don't care to do it.
%
% Ayn Rand

