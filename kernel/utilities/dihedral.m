% Computes the dihedral angle between vectors specified by
% the four sets of atomic coordinates. The atoms are assu-
% med to be bonded as A-B-C-D. Syntax:
%
%                 phi=dihedral(A,B,C,D)
%
% Parameters:
%
%     A  -  row vector of cartesian coordinates 
%           for atom A
%
%     B  -  row vector of cartesian coordinates 
%           for atom B
%
%     C  -  row vector of cartesian coordinates 
%           for atom C
%
%     D  -  row vector of cartesian coordinates 
%           for atom D
%
% Outputs:
%
%     phi - dihedral angle, degrees
%
% ilya.kuprov@weizmann.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=dihedral.m>

function phi=dihedral(A,B,C,D)

% Check consistency
grumble(A,B,C,D);

% Do the math
b1=(B-A)/norm(B-A,2); 
b2=(C-B)/norm(C-B,2);
b3=(D-C)/norm(D-C,2); 
phi=180*atan2(dot(norm(b2,2)*b1,cross(b2,b3)),...
              dot(cross(b1,b2),cross(b2,b3)))/pi;

end

% Consistency enforcement
function grumble(A,B,C,D)
if (~isnumeric(A))||(~isnumeric(B))||(~isnumeric(C))||(~isnumeric(D))||...
   (~isreal(A))||(~isreal(B))||(~isreal(C))||(~isreal(D))||...
   (numel(A)~=3)||(numel(B)~=3)||(numel(C)~=3)||(numel(D)~=3)||...
   (~isrow(A))||(~isrow(B))||(~isrow(C))||(~isrow(D))
    error('the arguments must be 3-element row vectors of real numbers.');
end
end

% Why is the plight of the poor not a matter of more sustained public
% discussion and more decisive government policy? [...] For many Ameri-
% cans, the real problem is not poverty at all; the real problem is the
% poor: [...] bad genes, bad work habits, and inadequate skills. Pover-
% ty is a symptom, a regrettable by-product of individual failings. The
% hardships experienced by the poor stem from their own shortcomings,
% not from any dysfunctions of the system; thus grand schemes to alle-
% viate poverty are inherently misguided. It might be appropriate, [...]
% for government to lend a modest helping hand, [...] but in the end,
% self-improvement, not social reform, is the only credible remedy.
% 
% Edward Royce

