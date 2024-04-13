% 2nd and 4th order Iserles product quadrature generators
% for one time propagation step in the case of state-inde-
% pendent Hamiltonian. Syntax:
%
%                   H=isergen(HL,HM,HR,dt)
%
% Parameters:
%
%     HL - Hamiltonian at the left edge of the interval
%
%     HM - [optional] Hamiltonian at the interval mid-
%          point; if this is empty, second order quad-
%          rature is used.
%
%     HR - Hamiltonian at the right edge of the interval
%
%     dt - interval duration, seconds
%
% Outputs:
%
%     H  - effective evolution generator, to be used
%          as exp(-1i*H*dt)
%
% i.kuprov@soton.ac.uk
% a.graham@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=isergen.m>

function H=isergen(HL,HM,HR,dt)

% Check consistency
grumble(HL,HM,HR,dt);

% Decide the quadrature
if isempty(HM)
    
    % Second order product quadrature
    H=(HL+HR)/2+(1i*dt/6)*(HL*HR-HR*HL);
    
else
    
    % Fourth order product quadrature
    H=(HL+4*HM+HR)/6+(1i*dt/12)*(HL*HR-HR*HL);
    
end

end

% Consistency enforcement
function grumble(HL,HM,HR,dt)
if (~isnumeric(HL))||(size(HL,1)~=size(HL,2))
    error('HL must be a square matrix');
end
if (~isnumeric(HR))||(size(HR,1)~=size(HR,2))
    error('HR must be a square matrix');
end
if ~isempty(HM)
    if (~isnumeric(HM))||(size(HM,1)~=size(HM,2))
        error('HM, if non-empty, must be a square matrix');
    end
end
if (~isnumeric(dt))||(~isreal(dt))||(~isscalar(dt))
    error('dt must be a real scalar');
end
end

% The hypothalamus plays a major role in the regulation of
% basic drives related to survival, including the so-called
% "four Fs": fighting, fleeing, feeding, and mating.
%
% Wayne Weiten, "Psychology: Themes and Variations"

