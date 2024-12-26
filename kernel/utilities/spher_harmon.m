% Spherical harmonics. Syntax:
%
%               Y=spher_harmon(l,m,theta,phi)
%
% Parameters:
%
%     l        - L quantum number
%
%     m        - M quantum number
%
%     theta    - an array of theta angles in radians
%
%     phi      - an array of phi angles in radians
%
% Outputs:
%
%     Y        - an array of spherical harmonics 
%                evaluated at the angles specified 
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=spher_harmon.m>

function Y=spher_harmon(l,m,theta,phi)

% Check consistency
grumble(l,m,theta,phi);

% Get Schmidt-normalized Legendres
S=legendre(l,cos(theta),'sch');
S=reshape(S,[l+1 numel(theta)]);
S=squeeze(S(abs(m)+1,:));
S=reshape(S,size(theta));

% Make spherical harmonics
if m==0
    Y=sqrt((2*l+1)/(4*pi))*S.*exp(1i*m*phi);
else
    Y=sqrt((2*l+1)/(4*pi))*S.*exp(1i*m*phi)/sqrt(2);
end

% Flip the sign if needed
if (m>0)&&(mod(m,2)==1), Y=-Y; end

end

% Consistency enforcement
function grumble(l,m,theta,phi)
if (~isnumeric(l))||(~isreal(l))||(~isscalar(l))||(mod(l,1)~=0)||(l<0)
    error('l must be a non-negative real integer.');
end
if (~isnumeric(m))||(~isreal(m))||(~isscalar(m))||(mod(m,1)~=0)||(m<-l)||(m>l)
    error('m must be a real integer from [-l,l] interval.');
end
if (~isnumeric(theta))||(~isreal(theta))
    error('theta must be numeric and real.');
end
if (~isnumeric(phi))||(~isreal(phi))
    error('phi must be numeric and real.');
end
end

% Don't pay any attention to what they write about 
% you. Just measure it in inches.
%
% Andy Warhol

