% Converts span and skew representation of a 3x3 interaction tensor
% (Herzfeld-Berger convention) into the corresponding matrix. Euler
% angles should be specified in radians. Syntax:
%
%                  M=spsk2mat(iso,sp,sk,alp,bet,gam)
%
% Parameters:
%
%        iso  - isotropic part of the interaction, defined as
%               (xx+yy+zz)/3 in terms of eigenvaues
%
%         sp  - interaction span, defined as the difference 
%               between the largest and the smallest eigenvalue
%
%         sk  - interaction skew, defined as 3*(yy-iso)/sp
%               where yy is the middle eigenvalue
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
% Note: the reverse transformation is ill-defined.
%
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=spsk2mat.m>

function M=spsk2mat(iso,sp,sk,alp,bet,gam)

% Check consistency
grumble(iso,sp,sk,alp,bet,gam);

% Compute eigenvalues
xx=iso-(1/6)*(3+sk)*sp;
yy=iso+sk*sp/3;
zz=iso+(1/6)*(3-sk)*sp;

% Rotate the matrix
R=euler2dcm(alp,bet,gam);
M=R*diag([xx yy zz])*R';

end

% Consistency enforcement
function grumble(iso,sp,sk,alp,bet,gam)
if (~isnumeric(iso))||(~isreal(iso))||(~isscalar(iso))||...
   (~isnumeric(sp))||(~isreal(sp))||(~isscalar(sp))||...
   (~isnumeric(sk))||(~isreal(sk))||(~isscalar(sk))||...
   (~isnumeric(alp))||(~isreal(alp))||(~isscalar(alp))||...
   (~isnumeric(bet))||(~isreal(bet))||(~isscalar(bet))||...
   (~isnumeric(gam))||(~isreal(gam))||(~isscalar(gam))
    error('all inputs must be real scalars.');
end
if abs(sk)>1
    error('skew cannot be outside [-1,+1] interval.');
end
if sp<0
    error('span must be positive.');
end
end

% Ilya,
% 
% you are being reprimanded for being a complete arsehole. However,
% I do not make this reprimand in public, despite its obviousness.
%
% Yours,
% Malcolm.
%
% Malcolm Levitt's email to IK, 08 Oct 2012.

