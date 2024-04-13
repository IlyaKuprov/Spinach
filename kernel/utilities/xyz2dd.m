% Converts coordinate specification of the dipolar interaction
% into the dipolar interaction constant, three Euler angles,
% and the dipolar interaction matrix. Syntax:
%
%      [d,alp,bet,gam,m]=xyz2dd(r1,r2,isotope1,isotope2)
%
% Parameters:
%
%    r1,2       - 3-element vectors of spin coordinates
%                 in Angstroms 
%
%    isotope1,2 - isotope specification strings, e.g.
%                 '13C'.
%
% Outputs:
%
%      d     - dipolar coupling constant, rad/s
%
%    alp     - alpha Euler angle, radians
%
%    bet     - beta Euler angle, radians
%
%    gam     - gamma Euler angle, radians
%
%      M     - dipolar interaction tensor, rad/s
%
% N.B. Euler angles are not uniquely defined for the orientati-
%      on of axial interactions (gamma angle can be anything).
%
% N.B. free-particle magnetogyric ratios are used, use xyz2hfc.m
%      if you have electrons in the system.
%
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=xyz2dd.m>

function [d,alp,bet,gam,M]=xyz2dd(r1,r2,isotope1,isotope2)

% Check consistency
grumble(r1,r2,isotope1,isotope2);

% Fundamental constants
hbar=1.054571628e-34; mu0=4*pi*1e-7;

% Get the distance
distance=norm(r2-r1,2);
        
% Get the ort
ort=(r2-r1)/distance;
        
% Get the dipolar interaction constant
d=spin(isotope1)*spin(isotope2)*hbar*mu0/(4*pi*(distance*1e-10)^3);
        
% Get the Euler angles
[alp,bet,~]=cart2sph(ort(1),ort(2),ort(3)); bet=pi/2-bet; gam=0;

% Compute the matrix if needed
if nargout>4

    % Get the dipolar coupling matrix
    M=d*[1-3*ort(1)*ort(1)   -3*ort(1)*ort(2)   -3*ort(1)*ort(3);
          -3*ort(2)*ort(1)  1-3*ort(2)*ort(2)   -3*ort(2)*ort(3);
          -3*ort(3)*ort(1)   -3*ort(3)*ort(2)  1-3*ort(3)*ort(3)];
  
    % Clean up rounding errors
    M=(M+M')/2; M=M-eye(3)*trace(M)/3;

end

end

% Consistency enforcement
function grumble(r1,r2,isotope1,isotope2)
if (~isnumeric(r1))||(~isreal(r1))||(numel(r1)~=3)
    error('r1 must be a three-element real vector.');
end
if (~isnumeric(r2))||(~isreal(r2))||(numel(r2)~=3)
    error('r2 must be a three-element real vector.');
end
if ~all(size(r1)==size(r2))
    error('r1 and r2 must have the same dimension.');
end
if (~ischar(isotope1))||(~ischar(isotope2))
    error('isotope specifications mst be character strings.');
end
end

% This principle is not a theorem, but a physical proposition, 
% that is, a vaguely stated and, strictly speaking, false as-
% sertion. Such assertions often happen to be fruitful sourc-
% es for mathematical theorems.
%
% Vladimir Arnold

