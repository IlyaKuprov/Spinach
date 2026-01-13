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
%               bottom up
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

function [states,coeffs]=enlev2ist(mult,lvl_num)

% Check consistency
grumble(mult,lvl_num);

% Energy level projector
P=zeros(mult,mult); 
P(mult-lvl_num+1,mult-lvl_num+1)=1;

% Spherical tensors
IST=irr_sph_ten(mult);

% Get all expansion coefficients and list all states
coeffs=cellfun(@(A)trace(A'*P)/trace(A'*A),IST);
states=transpose(0:(numel(IST)-1));

% Drop negligible states
idx=(abs(coeffs)>10*eps('double'));
coeffs=coeffs(idx); states=states(idx);

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

