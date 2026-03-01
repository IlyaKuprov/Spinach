% Exports phase-modulated optimal control waveforms into the
% format expected by Graham Smith's HiPER instrument. Syntax:
%
%           spinach2hiper(file_name,amp,phi,off,dt)
%
% Parameters:
%
%     file_name  - CSV file name, a character string
%                  without the extension
%
%     amp        - a vector of amplitudes
%
%     phi        - a vector of phases in radians
%
%     off        - transmitter offset in Hz
%
%     dt         - waveform slice duration, seconds
%
% Outputs:
%
%     this function writes a file
%
% Note: in practice, try both positive and negative pha-
%       ses - some instruments count phases clockwise, 
%       others counterclockwise.
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=spinach2hiper.m>

function spinach2hiper(file_name,amp,phi,off,dt)

% Check consistency
grumble(file_name,amp,phi,off,dt);

% Convert and wrap phases
phi=wrapTo360(180*phi(:)/pi);       % degrees [0 360]

% Slice timing and offsets
off=1e-6*repmat(off,size(phi)); % MHz
dt=1e9*repmat(dt,size(phi));    % ns
times=cumsum(dt)-dt;            % ns

% Make a table for export
pulse_table=table(times,off,phi,amp(:),'VariableNames',...
            {'time_ns','freq_MHz','phase_deg','amplitude'});

% Export the table into a CSV file
writetable(pulse_table,[file_name '.csv']);

end

% Consistency enforcement
function grumble(file_name,amp,phi,off,dt)
if ~ischar(file_name)
    error('file_name must be a character string.');
end
if (~isnumeric(amp))||(~isvector(amp))||(~isreal(amp))||any(~isfinite(amp),'all')
    error('amp must be a finite real vector.');
end
if (~isnumeric(phi))||(~isvector(phi))||(~isreal(phi))||any(~isfinite(phi),'all')
    error('phi must be a finite real vector.');
end
if numel(amp)~=numel(phi)
    error('amp and phi must have the same number of elements.');
end
if (~isnumeric(off))||(~isscalar(off))||(~isreal(off))||(~isfinite(off))
    error('off must be a finite real scalar.');
end
if (~isnumeric(dt))||(~isscalar(dt))||(~isreal(dt))||(~isfinite(dt))||(dt<=0)
    error('dt must be a positive finite real scalar.');
end
end

% It will get past reviewers because it invokes the PRL 
% Principle: make the paper as obscure as possible, then 
% ego prevents people from saying "I can't understand it"
% so it sneaks through. The only question is how high 
% the journal rank is going to be.
%
% Warren S. Warren, about a particularly egregious
% "NMR is Quantum Computing" ArXiV preprint.

