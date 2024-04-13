% Converts Euler angles (ZYZ active convention) into a direction
% cosine matrix. Syntax:
%
%                     R=euler2dcm(alpha,beta,gamma)
%
%                                 OR
%
%                    R=euler2dcm([alpha beta gamma])
%
% Parameters:
%
%    alpha,beta,gamma   -   Euler angles in radians (ZYZ
%                           active convention)
%
% Outputs:
%
%    R                  -   direction cosine matrix
%
% Note: the resulting rotation matrix is to be used as follows:
%
%         v=R*v    (for 3x1 vectors)
%         A=R*A*R' (for 3x3 interaction tensors)
%
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=euler2dcm.m>

function R=euler2dcm(arg1,arg2,arg3)

% Adapt to the input style
if nargin==1
    % Assume that a single input is a 3-vector
    alp=arg1(1); bet=arg1(2); gam=arg1(3);
elseif nargin==3
    % Assume that three inputs are Euler angles
    alp=arg1; bet=arg2; gam=arg3;
else
    % Bomb out in all other cases
    error('incorrect number of input arguments.');
end

% Check consistency
grumble(alp,bet,gam);

% Build the individual rotation matrices,
% as per Brink & Satchler, Fig 1a
R_alpha=[cos(alp)   -sin(alp)   0;          % (CCW around Z)
         sin(alp)    cos(alp)   0;
         0           0          1];
R_beta= [cos(bet)    0          sin(bet);   % (CCW around Y)
         0           1          0;
        -sin(bet)    0          cos(bet)];
R_gamma=[cos(gam)   -sin(gam)   0;          % (CCW around Z)
         sin(gam)    cos(gam)   0;
         0           0          1];
     
% Build the direction cosine matrix
R=R_alpha*R_beta*R_gamma;

end

% Consistency enforcement
function grumble(alp,bet,gam)
if (~isnumeric(alp))||(~isnumeric(bet))||(~isnumeric(gam))
    error('all inputs must be numeric.');
end
if (~isreal(alp))||(~isreal(bet))||(~isreal(gam))
    error('all inputs must be real.');
end
if (numel(alp)~=1)||(numel(bet)~=1)||(numel(gam)~=1)
    error('all inputs must have one element.');
end
end

% I am so hungry for any sight of anyone who's able to do 
% whatever it is he's doing!
%
% Ayn Rand, "Atlas Shrugged"

