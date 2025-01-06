% Returns a triangular waveform. Syntax:
%
%      waveform=triwave(amplitude,frequency,time_grid)
%
% Arguments:
%
%      amplitude  -  amplitude at the tooth top
%
%      frequency  -  waveform frequency, Hz
%
%      time_grid  -  vector of time points, seconds
%
% Outputs:
%
%      waveform   -  waveform array of the same shape
%                    as the time_grid input
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=triwave.m>

function waveform=triwave(amplitude,frequency,time_grid)

% Check consistency
grumble(amplitude,frequency,time_grid);

% Compute the waveform
waveform=abs(sawtooth(amplitude,frequency,time_grid));

end

% Consistency enforcement
function grumble(amplitude,frequency,time_grid)
if (numel(amplitude)~=1)||(~isnumeric(amplitude))||(~isreal(amplitude))
    error('amplitude parameter must be a real number.');
end
if (numel(frequency)~=1)||(~isnumeric(frequency))||(~isreal(frequency))
    error('frequency parameter must be a real number.');
end
if (~isnumeric(time_grid))||(~isreal(time_grid))
    error('time_grid parameter must be a vector of real numbers.');
end
end

% I am here to determine whether what you had just done is simple
% incompetence or deliberate sabotage.
%
% Joseph Stalin, to his generals, in May 1941.

