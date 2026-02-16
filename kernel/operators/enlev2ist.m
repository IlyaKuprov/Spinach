% Irreducible spherical tensor expansion of specific Zeeman 
% energy level projectors. Syntax:
%
%           [states,coeffs]=enlev2ist(mult,lvl_num)
%
% Parameters:
%
%     mult    - multipicity of the spin in question, a
%               positive integer
%
%     lvl_num - energy level number, counting from the
%               bottom up for spins and from top down
%               for bosons
%
%    particle - particle type, 'S' for a spin and 'B'
%               for a boson
%
% Outputs:
%
%     states  - states, in the Spinach IST basis index-
%               ing, that contribute to the operator in
%               question; use lin2lm to convert to L,M
%               spherical tensor indices
%
%     coeffs  - coefficients with which the ISTs enter
%               the linear combination
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=enlev2ist.m>

function [states,coeffs]=enlev2ist(mult,lvl_num,particle)

% Check consistency
grumble(mult,lvl_num);

% Energy level projector
P=zeros(mult,mult);

% Energy level counting
switch particle

    case 'S'

        % From the bottom up
        P(mult-lvl_num+1,mult-lvl_num+1)=1;

    case 'B'

        % From top to bottom
        P(lvl_num,lvl_num)=1;

    otherwise

        % Complain and bomb out
        error('unknown particle type.');

end

% Spherical tensor expansion
[states,coeffs]=oper2ist(P);

end

% Consistency enforcement
function grumble(mult,lvl_num)
if (~isnumeric(mult))||(~isscalar(mult))||...
   (~isreal(mult))||(mult<1)
    error('mult must be a positive integer.');
end
if (~isnumeric(lvl_num))||(~isscalar(lvl_num))||...
   (~isreal(lvl_num))||(lvl_num<1)||(lvl_num>mult)
    error('the specified Zeeman energy level does not exist.');
end
end

% |------------------------------------------------------------------------------------------------------------------|
% |   William de Wycombe, 1261         |   Ezra Pound, 1916                 |   Archibald Campbell, 1925             |
% |------------------------------------------------------------------------------------------------------------------|
% |   Sumer is icumen in,              |   Winter is icumen in,             |   Plumber is icumen in,                |
% |   Lhude sing cuccu.                |   Lhude sing goddamm.              |   Bludie big tu-du.                    |
% |   Groweth sed and bloweth med,     |   Raineth drop and staineth slop,  |   Bloweth lampe and showeth dampe,     |
% |   And springth the wde nu!         |   And how the wind doth ramm!      |   And dripth the wde thru!             |
% |                                    |                                    |                                        |
% |   Sing cuccu.                      |   Sing goddamm.                    |   Big tu-du.                           |
% |                                    |                                    |                                        |
% |   Awe bleteth after lomb,          |   Skiddeth bus and sloppeth us,    |   Thaweth drain and runneth bath,      |
% |   Lhouth after calue cu.           |   An ague hath my ham.             |   Saw saweth and scrueth scru.         |
% |   Bulluc sterteth, bucke uerteth,  |   Freezeth river, turneth liver,   |   Bull-kuk squirteth, leake spurteth,  |
% |   Murie sing cuccu.                |   Damn you, sing goddamm.          |   Wurry springth anew.                 |
% |------------------------------------------------------------------------------------------------------------------|

