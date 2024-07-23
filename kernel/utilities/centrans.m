% Central transition operators of half-integer spins in the
% Pauli basis. Syntax:
%
%                    A=centrans(mult,type)
%
% Parameters:
%
%     mult   - multipicity of the spin in question, an
%              even positive integer
%
%     type   - operator type: 'z' for polarisation, '+'
%              for raising, '-' for lowering
%
% Outputs:
%
%        A   - central transition operator
%
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=centrans.m>

function A=centrans(mult,type)

% Check consistency
grumble(mult,type);

% Build CT operator
A=spalloc(mult,mult,2);
switch type

    case 'x'

        % Sx on central transition
        A(mult/2,mult/2+1)=0.5;
        A(mult/2+1,mult/2)=0.5;

    case 'y'

        % Sy on central transition
        A(mult/2,mult/2+1)=-0.5i;
        A(mult/2+1,mult/2)=+0.5i;

    case 'z'

        % Sz on central transition
        A(mult/2,mult/2)=0.5;
        A(mult/2+1,mult/2+1)=-0.5;

    case '+'

        % S+ on central transition
        A(mult/2,mult/2+1)=1;

    case '-'

        % S- on central transition
        A(mult/2+1,mult/2)=1;

    otherwise

        % Complain and bomb out
        error('unknown CT operator type.');

end

% Make complex
A=complex(A);

end

% Consistency enforcement
function grumble(mult,type)
if (~isnumeric(mult))||(~isscalar(mult))||...
   (~isreal(mult))||(mult<2)||(mod(mult,2)~=0)
    error('mult must be an even positive integer.');
end
if (~ischar(type))||(~ismember(type,{'x','y','z','+','-'}))
    error('type must be ''x'', ''y'', ''z'', ''+'', or ''-''.');
end
end

% Prohibition agent Izzy Einstein bragged that he could find 
% liquor in any city under 30 minutes. In Chicago it took him
% 21 minutes, in Atlanta 17, and Pittsburgh just 11. But New
% Orleans set the record: 35 seconds. Einstein asked his taxi
% driver where to get a drink, and the driver handed him one.

