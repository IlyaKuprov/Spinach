% Converts fractional crystallographic coordinates to Cartesian
% coordinates. Syntax:
%
%      [XYZ,va,vb,vc]=frac2cart(a,b,c,alpha,beta,gamma,ABC)
%
% Parameters:
%
%    a,b,c            - three unit cell dimensions
%
%    alp,bet,gam      - three unit cell angles, degrees
%
%    ABC              - fractional atomic coordinates as
%                       Nx3 array of numbers
%
% Outputs:
%
%    XYZ              - Cartesian atomic coordinates as
%                       Nx3 array of numbers
%
%    va, vb, vc       - primitive lattice vectors
%
% ilya.kuprov@weizmann.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=frac2cart.m>

function [XYZ,va,vb,vc]=frac2cart(a,b,c,alp,bet,gam,ABC)

% Check consistency
grumble(a,b,c,alp,bet,gam,ABC);

% Compute the transformation matrix
v=a*b*c*sqrt(1-cosd(alp)^2-cosd(bet)^2-cosd(gam)^2+2*cosd(alp)*cosd(bet)*cosd(gam));
T=[a  b*cosd(gam)  c*cosd(bet);
   0  b*sind(gam)  c*(cosd(alp)-cosd(bet)*cosd(gam))/sind(gam);
   0  0            v/(a*b*sind(gam))];

% Apply the transformation matrix
XYZ=(T*ABC')';

% Get the primitive vectors
va=T(:,1); vb=T(:,2); vc=T(:,3);

end

% Consistency enforcement
function grumble(a,b,c,alp,bet,gam,ABC)
if (~isnumeric(a))||(~isscalar(a))||(~isreal(a))||(a<=0)||...
   (~isnumeric(b))||(~isscalar(b))||(~isreal(b))||(b<=0)||...
   (~isnumeric(c))||(~isscalar(c))||(~isreal(c))||(c<=0)
    error('a, b, c must be positive real numbers.');
end
if (~isnumeric(alp))||(~isscalar(alp))||(~isreal(alp))||...
   (~isnumeric(bet))||(~isscalar(bet))||(~isreal(bet))||...
   (~isnumeric(gam))||(~isscalar(gam))||(~isreal(gam))
    error('alp, bet, gam must be real scalars.');
end
if (~isnumeric(ABC))||(~isreal(ABC))||(size(ABC,2)~=3)
    error('ABC must be an Nx3 array of real numbers.');
end
end

% Socialism never took root in America because the poor see
% themselves not as an exploited proletariat but as tempora-
% rily embarrassed millionaires.
% 
% John Steinbeck

