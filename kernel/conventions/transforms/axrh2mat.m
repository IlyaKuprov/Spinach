% Converts axiality and rhombicity representation of a 3x3 interaction
% tensor into the corresponding matrix. Syntax:
%
%                 M=axrh2mat(iso,ax,rh,alp,bet,gam)
%
% Parameters:
%
%        iso  - isotropic part of the interaction, defined as
%               (xx+yy+zz)/3 in terms of eigenvaues
%
%         ax  - interaction axiality, defined as 2*zz-(xx+yy)
%               in terms of eigenvalues (Mehring order)
%
%         rh  - interaction rhombicity, defined as (yy-xx) in
%               terms of eigenvalues (Mehring order)
%
%        alp  - alpha Euler angle in radians
%
%        bet  - beta Euler angle in radians
%
%        gam  - gamma Euler angle in radians
%
% Outputs:
%
%          M  - 3x3 matrix
%
% Note: the inverse transformation is ill-defined.
%
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=axrh2mat.m>

function M=axrh2mat(iso,ax,rh,alp,bet,gam)

% Check consistency
grumble(iso,ax,rh,alp,bet,gam);

% Compute eigenvalues
xx=iso-(ax+3*rh)/6;
yy=iso-(ax-3*rh)/6;
zz=iso+ax/3;

% Rotate the molecule
R=euler2dcm(alp,bet,gam);
M=R*diag([xx yy zz])*R';

% Tidy up
M=(M+M')/2;

end

% Consistency enforcement
function grumble(iso,ax,rh,alp,bet,gam)
if (~isnumeric(iso))||(~isreal(iso))||(~isscalar(iso))||...
   (~isnumeric(ax))||(~isreal(ax))||(~isscalar(ax))||...
   (~isnumeric(rh))||(~isreal(rh))||(~isscalar(rh))||...
   (~isnumeric(alp))||(~isreal(alp))||(~isscalar(alp))||...
   (~isnumeric(bet))||(~isreal(bet))||(~isscalar(bet))||...
   (~isnumeric(gam))||(~isreal(gam))||(~isscalar(gam))
    error('all inputs must be real scalars.');
end
end

% I'm working to improve my methods, and every hour
% I save is an hour added to my life.
%
% Ayn Rand, "Atlas Shrugged"

