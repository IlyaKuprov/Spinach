% Kronecker symbol. Syntax:
%
%                      d=krondelta(a,b)
%
% Parameters:
%
%       a  - an integer number
%
%       b  - an integer number
%
% Output:
%
%       d  - a logical number
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=krondelta.m>

function d=krondelta(a,b)

% Check consistency
grumble(a,b);

% Compute the answer
if a==b, d=true(); else, d=false(); end

end

% Consistency enforcement
function grumble(a,b)
if (~isnumeric(a))||(~isnumeric(b))||(~isscalar(a))||(~isscalar(b))||...
   (~isreal(a))||(~isreal(b))||(mod(a,1)~=0)||(mod(b,1)~=0)
    error('a and b must be real integer scalars.');
end
end

% Die ganzen Zahlen hat der liebe Gott gemacht,
% alles andere ist Menschenwerk.
%
% Leopold Kronecker

