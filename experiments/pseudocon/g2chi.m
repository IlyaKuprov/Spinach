% Calculates a high-termperature estimate of the magnetic suscep-
% tibility tensor from the user-supplied g-tensor. Syntax:
%
%                          chi=g2chi(g,T,S)
%
% Parameters:
%
%      g - 3x3 g-tensor matrix in Bohr magneton units
%
%      T - absolute temperature in Kelvin
%
%      S - electron spin (1/2, 1, 3/2, etc.)
%
% Outputs:
%
%    chi - 3x3 magnetic susceptibility tensor in cubic Angstrom
%
% ilya.kuprov@weizmann.ac.uk
% e.suturina@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=g2chi.m>

function chi=g2chi(g,T,S)

% Check consistency
grumble(g,T,S);

% Fundamental constants
mu_b=9.274009994e-24;
mu_0=4*pi*1e-7;
k_b=1.38064852e-23;

% Diagonalise the g-tensor
[V,D]=eig(g);

% Apply the Curie law to the eigenvalues
D=S*(S+1)*mu_0*(mu_b^2)*(D.^2)/(3*k_b*T);

% Compose the susceptibility tensor
chi=1e30*V*D*inv(V); %#ok<MINV>

end

% Consistency enforcement
function grumble(g,T,S)
if (~isnumeric(g))||(~isreal(g))||...
   (~ismatrix(g))||any(size(g)~=[3 3])
    error('g must be a real 3x3 matrix.');
end
if (~isnumeric(T))||(~isreal(T))||...
   (~isscalar(T))||(T<=0)
    error('T must be a positive real scalar.');
end
if (~isnumeric(S))||(~isreal(S))||...
   (~isscalar(S))||(S<=0)||(mod(2*S,1)~=0)
    error('S must be a positive integer or half-integer number.');
end
end

% Sometimes it pays to stay in bed on Monday, rather than
% spending the rest of the week debugging Monday's code.
%
% Christopher Thompson

