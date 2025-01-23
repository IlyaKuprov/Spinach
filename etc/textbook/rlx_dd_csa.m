% Redfield theory expressions for some relaxation and cross-
% relaxation rates in a CSA-DD-CSA system with two spin-1/2
% particles. Syntax:
%
%    [A,B,X]=rlx_dd_csa(B0,tau_c,isotopes,deltas,coords)
%
% Parameters:
%
%    B0       - magnet field, Tesla
%
%    tau_c    - rotational correlation time, seconds
%
%    isotopes - the spins involved, e.g. {'13C','19F'}
%
%    deltas   - a cell array with two symmetric 3x3
%               chemical shift tensors in ppm
%
%    coords   - a cell array with two 1x3 Cartesi-
%               an coordinate vectors in Angstrom
%
% Outputs:
%
%    (A,B).r(1,2).csa   - CSA contribution to R1 and R2
%                         rates of spins A and B
%
%    (A,B).r(1,2).dd    - dipolar contribution to R1 and 
%                         R2 rates of spins A and B
%
%    (A,B).r(1,2).total - total R1 and R2 rates
%
%    (A,B).trosy.dd     - dipole contribution to the 
%                         transverse relaxation rate
%                         of TROSY doublet components
%                         of spins A and B
%
%    (A,B).trosy.csa    - CSA contribution to the 
%                         transverse relaxation rate
%                         of TROSY doublet components
%                         of spins A and B
%
%    (A,B).trosy.xc     - cross-correlation contribu-
%                         tion to the transverse rela-
%                         xation rate of TROSY doublet
%                         components of spins A and B
%
%    (A,B).trosy.total_bro - transverse relaxation rate of
%                            the broad TROSY doublet com-
%                            ponent of spins A and B
%
%    (A,B).trosy.total_nar - transverse relaxation rate of
%                            the narrow TROSY doublet com-
%                            ponent of spins A and B
%
%     X - longitudinal cross-relaxation rate between 
%         spins A and B
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=rlx_dd_csa.m>

function [A,B,X]=rlx_dd_csa(B0,tau_c,isotopes,deltas,coords)

% Check consistency
grumble(B0,tau_c,isotopes,deltas,coords);

% Carrier frequencies
omegaA=B0*(1+trace(1e-6*deltas{1})/3)*spin(isotopes{1});
omegaB=B0*(1+trace(1e-6*deltas{2})/3)*spin(isotopes{2});

% Spectral power density function
J=@(omega)tau_c/(1+omega^2*tau_c^2);

% Aniso parts of Zeeman tensors
ZA=1e-6*B0*(deltas{1}-eye(3)*trace(deltas{1})/3)*spin(isotopes{1});
ZB=1e-6*B0*(deltas{2}-eye(3)*trace(deltas{2})/3)*spin(isotopes{2});

% Dipole-dipole coupling tensor
[~,~,~,~,DD]=xyz2dd(coords{1},coords{2},isotopes{1},isotopes{2});

% Blicharski invariants and products
[~,DsqZA]=blinv(ZA); [~,DsqZB]=blinv(ZB); [~,DsqDD]=blinv(DD);
[~,X_DD_ZA]=blprod(DD,ZA); [~,X_DD_ZB]=blprod(DD,ZB);

% Longitudinal relaxation rates
A.r1.csa=(2/15)*DsqZA*J(omegaA);
A.r1.dd=(1/90)*DsqDD*(3*J(omegaA)+J(omegaA-omegaB)+6*J(omegaA+omegaB));
A.r1.total=A.r1.csa+A.r1.dd;
B.r1.csa=(2/15)*DsqZB*J(omegaB);
B.r1.dd=(1/90)*DsqDD*(3*J(omegaB)+J(omegaB-omegaA)+6*J(omegaB+omegaA));
B.r1.total=B.r1.csa+B.r1.dd;

% Transverse relaxation rates
A.r2.csa=(1/45)*DsqZA*(4*J(0)+3*J(omegaA));
A.r2.dd=(1/180)*DsqDD*(4*J(0)+3*J(omegaA)+J(omegaA-omegaB)+...
                       6*J(omegaB)+6*J(omegaA+omegaB));
A.r2.total=A.r2.csa+A.r2.dd;
B.r2.csa=(1/45)*DsqZB*(4*J(0)+3*J(omegaB));
B.r2.dd=(1/180)*DsqDD*(4*J(0)+3*J(omegaB)+J(omegaB-omegaA)+...
                       6*J(omegaA)+6*J(omegaB+omegaA));
B.r2.total=B.r2.csa+B.r2.dd;

% TROSY component relaxation rates
A.trosy.dd=(1/180)*DsqDD*(4*J(0)+3*J(omegaA)+J(omegaA-omegaB)+...
                          3*J(omegaB)+6*J(omegaA+omegaB));
A.trosy.csa=(1/45)*DsqZA*(4*J(0)+3*J(omegaA))+(1/15)*DsqZB*J(omegaB);                     
A.trosy.xc=(1/45)*X_DD_ZA*(4*J(0)+3*J(omegaA));
A.trosy.total_bro=A.trosy.dd+A.trosy.csa+abs(A.trosy.xc);
A.trosy.total_nar=A.trosy.dd+A.trosy.csa-abs(A.trosy.xc);

B.trosy.dd=(1/180)*DsqDD*(4*J(0)+3*J(omegaB)+J(omegaB-omegaA)+...
                          3*J(omegaA)+6*J(omegaB+omegaA));
B.trosy.csa=(1/45)*DsqZB*(4*J(0)+3*J(omegaB))+(1/15)*DsqZA*J(omegaA);                     
B.trosy.xc=(1/45)*X_DD_ZB*(4*J(0)+3*J(omegaB));
B.trosy.total_bro=B.trosy.dd+B.trosy.csa+abs(B.trosy.xc);
B.trosy.total_nar=B.trosy.dd+B.trosy.csa-abs(B.trosy.xc);

% Cross-relaxation rate
X=(1/90)*DsqDD*(6*J(omegaA+omegaB)-J(omegaA-omegaB));
        
end

% Consistency enforcement
function grumble(B0,tau_c,isotopes,deltas,coords)
if (~isnumeric(B0))||(~isreal(B0))||(~isscalar(B0))
    error('B0 must be a real number.');
end
if (~isnumeric(tau_c))||(~isreal(tau_c))||...
   (~isscalar(tau_c))||(tau_c<=0)
    error('tau_c must be a positive real number.');
end
if ~iscell(isotopes)
    error('isotopes must be a cell array of character strings.');
end
if ~iscell(deltas)
    error('deltas must be a cell array of 3x3 matrices.');
end
if (~issymmetric(deltas{1}))||(~issymmetric(deltas{2}))
    error('this function only supports symmetric shift tensors.');
end
if ~iscell(coords)
    error('coords must be a cell array of 1x3 vectors.');
end
[~,mult_a]=spin(isotopes{1}); [~,mult_b]=spin(isotopes{2});
if (mult_a~=2)||(mult_b~=2)
    error('this function only supports spin-1/2 particles.');
end
end

% "We show that LaTeX users were slower than Word users, 
%  wrote less text in the same amount of time, and pro-
%  duced more typesetting, orthographical, grammatical,
%  and formatting errors. On most measures, expert LaTeX
%  users performed even worse than novice Word users."
%
% https://doi.org/10.1371/journal.pone.0125830

