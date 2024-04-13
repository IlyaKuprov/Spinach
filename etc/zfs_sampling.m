% Gadolinium ZFS probability distribution function for DOTA-type 
% ligand complexes in cryogenic water-methanol glasses. The para-
% meters match those given in Figure 5 of 
%
%                https://doi.org/10.1007/BF03166762
%
% Syntax:
%
%          [D,E,W]=zfs_sampling(npoints_d,npoints_e,tol)
%
% Parameters:
%
%      npoints_d - number of Gauss-Legendre quadrature
%                  points in D
%
%      npoints_e - number of Gauss-Legendre quadrature 
%                  points in E
%
%      tol       - tolerance for integration weights 
%                  below which grid points are dropped
%
% Outputs:
%
%      D - a vector of D values at each integration
%          grid point
%
%      E - a vector of E values at each integration
%          grid point
%
%      W - a vector of weights for each integration
%          grid point
%
% Notes: the function also creates a figure with the distributi-
%        ons it has used for D and E parameters.
%
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=zfs_sampling.m>

function [D,E,W]=zfs_sampling(npoints_d,npoints_e,tol)

% Check consistency
grumble(npoints_d,npoints_e,tol);

% Generate Gauss-Legendre point set for D/D1
[X,WX]=gaussleg(-2,2,npoints_d);

% Get the standard deviation for unit FWHM
sigma=1/(2*sqrt(2*log(2)));

% Refract weights through a double Gaussian
WX=WX.*(normpdf(X,-1,sigma)+normpdf(X,+1,sigma)); WX=WX/sum(WX);

% Plot the double Gaussian
figure(); scale_figure([1.50 0.75]); subplot(1,2,1); 
plot(X,normpdf(X,-1,sigma)+normpdf(X,+1,sigma),'r-');
ktitle('$D/D_0$ distribution'); kxlabel('$D/D_0$'); 
kylabel('probability density'); xlim tight; kgrid;

% Generate Gauss-legendre point set for E/D
[Y,WY]=gaussleg(0,1/3,npoints_e);

% Refract weights through a quadratic function
WY=WY.*(-(Y-0.25).^2+0.0625); WY=WY/sum(WY);

% Plot the quadratic function
subplot(1,2,2); plot(Y,-(Y-0.25).^2+0.0625,'r-');
ktitle('$E/D$ distribution'); kxlabel('$E/D$'); 
kylabel('probability density'); xlim tight; kgrid;

% Kron the weights
D=kron(X,ones(size(Y)));
E=D.*kron(ones(size(X)),Y);
W=kron(WX,WY); W=W/sum(W);

% Ignore small weights
D(W<tol)=[]; E(W<tol)=[]; W(W<tol)=[];

end

% Consistency enforcement
function grumble(npoints_d,npoints_e,tol)
if (~isnumeric(npoints_d))||(~isreal(npoints_d))||...
   (npoints_d<5)||(mod(npoints_d,1)~=0)
    error('npoints_d must be a real integer greater than 5.');
end
if (~isnumeric(npoints_e))||(~isreal(npoints_e))||...
   (npoints_e<5)||(mod(npoints_e,1)~=0)
    error('npoints_e must be a real integer greater than 5.');
end
if (~isnumeric(tol))||(~isreal(tol))||(tol<0)
    error('tol must be a positive real number much smaller than 1.');
end
end

% It is a cliche that most cliches are true, but 
% then, like most cliches, that cliche is untrue.
%
% Stephen Fry

