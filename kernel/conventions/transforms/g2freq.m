% Converts g-tensor units into electron Zeeman frequency 
% units. Syntax:
%
%                        f=g2freq(g,B)
%
% Parameters:
%
%     g  -  frequency in g-tensor units (Bohr magnetons)
%
%     B  -  magnet field in Tesla
%
% Outputs:
%
%     f  -  frequency in Hz
%
% ilya.kuprov@weizmann.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=g2freq.m>

function f=g2freq(g,B)

% Check consistency
grumble(g,B)

% Get the free electron carrier frequency
omega=B*spin('E')/(2*pi);

% Scale the frequency
f=g.*omega/2.0023193043622;

end

% Consistency enforcement
function grumble(g,B)
if (~isnumeric(g))||(~isreal(g))||(~isnumeric(B))||(~isreal(B))
    error('both arguments must be numeric and real.');
end
end

% "Authors are listed in order of degree of belief in
%  the central thesis."
%
% A footnote in http://www.jstor.org/stable/3328150

