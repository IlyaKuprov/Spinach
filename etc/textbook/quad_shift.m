% Second order shift of the centre of gravity of the powder pattern
% of |S,m> to |S,m-1> transition in the NMR spectrum of a quadrupo-
% lar nucleus with spin S. Equation (3) from
%
%           https://doi.org/10.1016/0009-2614(85)85414-2
%
% Syntax:
%
%                   delta=quad_shift(Cq,eta,v0,S,m)
%
% Parameters:
%
%    Cq    - quadrupolar constant, Hz
%
%    eta   - quadrupolar asymmetry parameter
%
%    v0    - Larmor frequency of the nucleus, Hz
%
%    S     - spin quantum number of the nucleus
%
%    m     - projection quantum number of the
%            starting energy level
%
% Outputs:
%
%    delta - quadrupolar shift in ppm
%
% Note: a few papers contain an incorrect version of this expressi-
%       on; the one used here was tested against pure numerics and
%       found to be correct.
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=quad_shift.m>

function delta=quad_shift(Cq,eta,v0,S,m)

% Check consistency
grumble(Cq,eta,v0,S,m);

% Use Samoson's expression
delta=-1e6*(3/40)*(Cq/v0)^2*(1+eta^2/3)*...
           (S*(S+1)-9*m*(m-1)-3)/(S^2*(2*S-1)^2);

end

% Consistency enforcement
function grumble(Cq,eta,v0,S,m)
if (~isnumeric(Cq))||(~isreal(Cq))||(~isscalar(Cq))
    error('Cq must be a real scalar.');
end
if (~isnumeric(eta))||(~isreal(eta))||(~isscalar(eta))
    error('eta must be a real scalar.');
end
if (~isnumeric(v0))||(~isreal(v0))||(~isscalar(v0))
    error('v0 must be a real scalar.');
end
if (~isnumeric(S))||(~isreal(S))||(~isscalar(S))||...
   (S<1)||((mod(S,1)~=0)&&(mod(2*S,1)~=0))
    error('S must be integer or half-integer, and greater than 1/2.');
end
if (~isnumeric(m))||(~isreal(m))||(~isscalar(m))||...
   (m>S)||(m<(1-S))||((mod(m,1)~=0)&&(mod(2*m,1)~=0))
    error('m must indicate an existing |S,m> to |S,m-1> transition.');
end
end

% Two computer games have influenced IK in a profound 
% way: "Baldur's Gate 2" (by the story of Jon Irenicus)
% and "Planescape: Torment" (by many things, but most
% of all by the story of Deionarra). Do not start these
% games unless you have a month to spare.

