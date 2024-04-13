% Bosonic operators; ad hoc numbering scheme until such time
% as we code up Weyl algebras. Syntax:
%
%                   A=boson_oper(mult,opspec)
%
% Parameters:
%
%   mult    - number of population levels to consider in
%             the bosonic mode, an integer
%
%   opspec  - operator specification, an integer: -1 for
%             creation operator, -2 for population number
%             operator, -3 for annihilation operator, -4
%             for the empty mode operator
%
% Outputs:
%
%   A       - a [mult]x[mult] matrix
%
% Note: arrays are declared complex at creation to avoid 
%       expensive reallocation operations later on.
%
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=boson_oper.m>

function A=boson_oper(mult,opspec)

% Check consistency
grumble(mult,opspec);

% Decide the type
switch opspec

    case -1

        % Creation operator
        diags=sqrt(1:mult); A=spdiags(diags',-1,mult,mult);

    case -2

        % Number operator
        diags=0:(mult-1); A=spdiags(diags',0,mult,mult);

    case -3

        % Annihilation operator
        diags=sqrt(0:(mult-1)); A=spdiags(diags',+1,mult,mult);

    case -4

        % Empty cavity operator
        A=spalloc(mult,mult,1); A(1,1)=1;
        
    otherwise

        % Complain and bomb out
        error('unknown operator type.');

end

% Efficiency
A=complex(A);

end

% Consistency enforcement
function grumble(mult,opspec)
if (~isnumeric(mult))||(~isscalar(mult))||(~isreal(mult))||...
   (mult<1)||(mod(mult,1)~=0)
    error('mult must be a real and positive integer.');
end
if (~isnumeric(opspec))||(~isscalar(opspec))||...
   (~isreal(opspec))||(mod(opspec,1)~=0)
    error('opspec must be a real integer.');
end
end

% There are only two kinds of languages: the ones people
% complain about and the ones nobody uses.
%
% Bjarne Stroustrup,
% creator of C++

