% Wigner D matrices, defined as (Brink & Satchler, Eq 2.13):
%
%       D=expm(-1i*Lz*alp)*expm(-1i*Ly*bet)*expm(-1i*Lz*gam);
%
% where Lx, Ly, Lz are Pauli matrices. Syntax:
%
%                      D=wigner(l,alp,bet,gam)
%
% Parameters:
%
%       l    - rank of the Wigner matrix, may be half-integer
%
%     alp    - alp Euler angle, radians
%
%     bet    - bet Euler angle, radians
%
%     gam    - gam Euler angle, radians
%
% ZYZ convention is used for Euler angles, see Brink and Satchler,
% Figures 1 and 2. 
%
% Outputs: 
%
%       D    - Wigner D matrix with rows and columns sorted
%              by descending ranks, for example (l=2):
%    
%                      [D( 2,2)  ...  D( 2,-2)
%                         ...    ...    ...  
%                       D(-2,2)  ...  D(-2,-2)]
%
% The output is to be used as y=D*x, where x is a column vector of
% irreducible spherical tensor coefficients, listed vertically in
% the order: T(2,2), T(2,1), T(2,0), T(2,-1), T(2,-2). 
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=wigner.m>

function D=wigner(l,alp,bet,gam)

% Check consistency
grumble(l,alp,bet,gam);

% Rank 1 hard-coded for speed
if l==1
    
    % Blob components
    egp1=exp(-1i*1*gam); egm1=exp(-1i*(-1)*gam); 
    eap1=exp(-1i*1*alp); eam1=exp(-1i*(-1)*alp);
    sb=sin(bet); cb=cos(bet);
    
    % This ugly blob runs faster than the generic expm call below
    D=[eap1*egp1*(1/2+cb/2), -eap1*sb/sqrt(2), eap1*egm1*(1/2-cb/2); ...
       egp1*sb/sqrt(2), cb, -egm1*sb/sqrt(2); eam1*egp1*(1/2-cb/2),  ...
       eam1*sb/sqrt(2), eam1*egm1*(1/2+cb/2)];
    
% Rank 2 hard-coded for speed
elseif l==2
    
    % Blob components
    egp2=exp(-1i*2*gam); egp1=exp(-1i*1*gam); 
    egm1=exp(-1i*(-1)*gam); egm2=exp(-1i*(-2)*gam);
    eap2=exp(-1i*2*alp); eap1=exp(-1i*1*alp);
    eam1=exp(-1i*(-1)*alp); eam2=exp(-1i*(-2)*alp);
    sb=sin(bet); cb=cos(bet); cbh4=cos(bet/2)^4; sbh4=sin(bet/2)^4;
    
    % This ugly blob runs faster than the generic expm call below
    D=[eap2*cbh4*egp2, eap2*(-sb*(1+cb)/2)*egp1, eap2*sqrt(3/8)*sb^2,        ...
       eap2*(-sb*(1-cb)/2)*egm1, eap2*sbh4*egm2; eap1*(sb*(cb+1)/2)*egp2,    ...
       eap1*(2*cb-1)*(1+cb)/2*egp1, eap1*(-sqrt(3/2)*sb*cb),                 ...
       eap1*(2*cb+1)*(1-cb)/2*egm1, eap1*(sb*(cb-1)/2)*egm2;                 ...
       sqrt(3/8)*sb^2*egp2, sqrt(3/2)*sb*cb*egp1, (3*cb^2-1)/2,              ...
       -sqrt(3/2)*sb*cb*egm1, sqrt(3/8)*sb^2*egm2; -eam1*(sb*(cb-1)/2)*egp2, ...
       eam1*(2*cb+1)*(1-cb)/2*egp1, eam1*(sqrt(3/2)*sb*cb),                  ...
       eam1*(2*cb-1)*(1+cb)/2*egm1, eam1*(-sb*(cb+1)/2)*egm2;                ...
       eam2*sbh4*egp2, eam2*sb*(1-cb)/2*egp1, eam2*sqrt(3/8)*sb^2,           ...
       eam2*sb*(1+cb)/2*egm1, eam2*cbh4*egm2];
    
else
    
    % Get Pauli matrices
    L=pauli(2*l+1);
    
    % Compute Wigner matrix (Brink & Satchler, Eq 2.13)
    D=expm(-1i*L.z*alp)*expm(-1i*L.y*bet)*expm(-1i*L.z*gam);

end

end

% Consistency enforcement
function grumble(l,alp,bet,gam)
if (~isnumeric(l))||(~isreal(l))||(~isscalar(l))||(mod(l,1/2)~=0)||(l<0)
    error('l must be a non-negative real integer or half-integer.');
end
if (~isnumeric(alp))||(~isreal(alp))||(~isscalar(alp))||...
   (~isnumeric(bet))||(~isreal(bet))||(~isscalar(bet))||...
   (~isnumeric(gam))||(~isreal(gam))||(~isscalar(gam))
    error('alp, bet and gam must be real scalars.');
end
end

% Every stink that fights the ventilator thinks it
% is Don Quixote.
%
% Stanislaw Lec

