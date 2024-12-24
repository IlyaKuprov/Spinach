% Amplitude envelopes of pulse waveforms. Syntax:
%
%          waveform=pulse_shape(pulse_name,npoints)
%
% Parameters:
%
%     pulse_name - the name of the pulse (see function text)
%
%     npoints    - number of points in the pulse
%
% Output:
%
%     waveform - normalised waveform as a vector
%
% ilya.kuprov@weizmann.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=pulse_shape.m>

function waveform=pulse_shape(pulse_name,npoints)

% Check consistency
grumble(pulse_name,npoints);

% Choose the shape
switch pulse_name
    
    case 'gaussian'
        
        time_grid=linspace(-2,2,npoints);
        waveform=normpdf(time_grid)/sqrt(2);
        
    case 'sinc5'
        
        time_grid=linspace(-5,5,npoints);
        waveform=pi*sinc(time_grid);
        
    case 'sinc3'
        
        time_grid=linspace(-3,3,npoints);
        waveform=pi*sinc(time_grid);
        
    case 'rectangular'
        
        waveform=ones(1,npoints);
        
    otherwise
        
        % Complain and bomb out
        error('unknown pulse name.');
        
end

end

% Consistency enforcement
function grumble(pulse_name,npoints)
if ~ischar(pulse_name)
    error('pulse_name parameter must be a character string.');
end
if (numel(npoints)~=1)||(~isnumeric(npoints))||(~isreal(npoints))||...
   (npoints<1)||(mod(npoints,1)~=0)
    error('npoints parameter must be a positive real integer greater than 1.');
end
end

% Let me start with a parable. It concerns an Eastern European country
% whose parliament was considering a total smoking ban. In response, a
% consortium of tobacco companies demonstrated that the savings made in
% healthcare as a result of the decline in smoking-related diseases were
% chicken feed besides the reduced payout in pensions as the result of
% premature death - not to mention the fiscal increment from the habit.
%
% Stewart Dakers

