% Converts g-tensor units into electron Zeeman frequency 
% units. Syntax:
%
%                        f=g2freq(g,B)
%
% Parameters:
%
%     g  -  g-values, scalar or array
%
%     B  -  magnetic field in Tesla
%
% Outputs:
%
%     f  -  frequency in Hz
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=g2freq.m>

function f=g2freq(g,B)

% Check consistency
grumble(g,B);

% Get the free electron carrier frequency
omega=B*spin('E')/(2*pi);

% Scale the frequency
f=g.*omega/2.0023193043622;

end

% Consistency enforcement
function grumble(g,B)
if (~isnumeric(g))||(~isreal(g))||any(~isfinite(g),'all')
    error('g must be a finite real numeric array.');
end
if (~isnumeric(B))||(~isscalar(B))||(~isreal(B))||(~isfinite(B))
    error('B must be a finite real scalar.');
end
end

% "Authors are listed in order of degree of belief in
%  the central thesis."
%
% A footnote in http://www.jstor.org/stable/3328150

