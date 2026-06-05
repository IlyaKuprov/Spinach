% Instantaneous frequency trajectory from a complex time-domain
% signal by regularised phase differentiation. Syntax:
%
%          freq=inst_freq(signal,dt,npoints,poly_order)
%
% Parameters:
%
%    signal - complex row or column vector 
%             with time-domain signal
%
%    dt     - time step duration between 
%             signal points (seconds)
%
%    npoints    - odd number of signal points in the
%                 local least-squares window
%
%    poly_order - local polynomial order used for
%                 phase differentiation
%
% Outputs:
%
%    freq   - instantaneous frequency trajec-
%             tory (Hz), same size as signal
%             and on the same time grid
%
% Note: signal phase is unwrapped first and then differentiated
%       by Savitzky-Golay local least-squares polynomial fits.
%       This regularises numerical phase noise before the deri-
%       vative is taken.
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=inst_freq.m>

function freq=inst_freq(signal,dt,npoints,poly_order)

% Check consistency
grumble(signal,dt,npoints,poly_order);

% Input into column
signal_col=signal(:);

% Unwrap the phase trajectory
phase_col=unwrap(angle(signal_col));

% Differentiate the phase using local least-squares fits
freq=sgolaydiff(phase_col,1,npoints,poly_order)/(2*pi*dt);

% Shape back to input shape
freq=reshape(freq,size(signal));

end

% Consistency enforcement
function grumble(signal,dt,npoints,poly_order)
if (~isnumeric(signal))||(~isvector(signal))||isempty(signal)||isreal(signal)
    error('signal must be a non-empty complex numeric vector.');
end
if numel(signal)<3
    error('signal must contain at least three points.');
end
if any(~isfinite(signal))
    error('signal must not contain Inf or NaN.');
end
if (~isnumeric(dt))||(~isreal(dt))||(~isscalar(dt))||(~isfinite(dt))||(dt<=0)
    error('dt must be a positive real scalar.');
end
if (~isnumeric(npoints))||(~isreal(npoints))||(numel(npoints)~=1)||...
   (npoints<3)||(mod(npoints,1)~=0)
    error('npoints must be a positive odd integer greater than two.');
end
if mod(npoints,2)~=1
    error('npoints must be an odd integer.');
end
if npoints>numel(signal)
    error('npoints must not exceed the number of signal points.');
end
if (~isnumeric(poly_order))||(~isreal(poly_order))||(numel(poly_order)~=1)||...
   (poly_order<1)||(mod(poly_order,1)~=0)
    error('poly_order must be a positive real integer.');
end
if poly_order>=npoints
    error('poly_order must be smaller than npoints.');
end
end

% Sir Isaac Newton was rigorously puritanical: when one of
% his few friends told him "a loose story about a nun", he
% ended their friendship. He is not known to have ever had
% a romantic relationship of any kind, and is believed to
% have died a virgin.

