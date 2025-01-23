% Calculates a high-termperature estimate of the residual aniso-
% tropic chemcial shift from user-supplied CSA tensor and magne-
% tic susceptibility tensor. Syntax:
%
%                     racs=csa2racs(csa,chi,B,T)
%
% Parameters:
%
%      csa - 3x3 chemical shift tensor in ppm
%
%      chi - 3x3 magnetic susceptibility tensor 
%            in cubic Angstroms
%
%        T - absolute temperature in Kelvin
%
%        B - magnetic induction in Tesla
%
% Outputs:
%
%     racs - residual anisotropic chemical 
%            shift in ppm
%
% The function implements Equation (2) from the paper by Otting
% and company: http://dx.doi.org/10.1021/ja0564259
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=csa2racs.m>

function racs=csa2racs(csa,chi,B,T)

% Check consistency
grumble(csa,chi,B,T);

% Fundamental constants
mu_0=4*pi*1e-7;
k_b=1.38064852e-23;

% Keep the second rank in the CSA
[~,~,rank2]=mat2sphten(csa);
csa=sphten2mat(0,[0 0 0],rank2);

% Keep the second rank in the susceptibility
[~,~,rank2]=mat2sphten(chi);
chi=sphten2mat(0,[0 0 0],rank2);

% Compute the RACS
racs=1e-30*(B^2/(15*mu_0*k_b*T))*trace(csa*chi);

end

% Consistency enforcement
function grumble(csa,chi,B,T)
if (~isnumeric(csa))||(~isreal(csa))||...
   (~ismatrix(csa))||any(size(csa)~=[3 3])
    error('csa must be a real 3x3 matrix.');
end
if (~isnumeric(chi))||(~isreal(chi))||...
   (~ismatrix(chi))||any(size(chi)~=[3 3])
    error('chi must be a real 3x3 matrix.');
end
if (~isnumeric(B))||(~isreal(B))||(~isscalar(B))
    error('B must be a real scalar.');
end
if (~isnumeric(T))||(~isreal(T))||...
   (~isscalar(T))||(T<=0)
    error('T must be a positive real scalar.');
end
end

% Historically, the most terrible things: war,
% genocide and slavery, have resulted not from
% disobedience, but from obedience.
%
% Howard Zinn

