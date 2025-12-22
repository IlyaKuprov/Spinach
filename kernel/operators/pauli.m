% Pauli spin operators (sparse, see below for normalisa-
% tion conventions) for a spin of a user-specified ener-
% gy level multiplicity. Syntax:
%
%                      S=pauli(mult)
%
% Parameters:
%
%     mult - an integer specifying the 
%            multiplicity of the spin
%
% Outputs:
%
%     S.u - unit operator
%
%     S.p - raising operator
%
%     S.m - lowering operator
%
%     S.x - Sx observable operator
%
%     S.y - Sy observable operator
%
%     S.z - Sz observable operator
%
% Note: the matrices are normalised to obey the following
%       commutation relations for all multiplicities: 
%
%                      [S.x,S.y]=1i*S.z
%                      [S.y,S.z]=1i*S.x
%                      [S.z,S.x]=1i*S.y
%
% Note: raising and lowering operators are defined as:
%
%                      S.p=S.x+1i*S.y
%                      S.m=S.x-1i*S.y
%
% Note: arrays are declared complex at creation to avoid 
%       expensive reallocation operations later on.
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=pauli.m>

function S=pauli(mult)

% Ensure internal consistency
grumble(mult); mult=double(mult);

% Get Pauli matrices
if mult==2
    
    % Spin-half matrices are hard-coded for speed
    S.u=sparse(complex([1 0; 0 1]));
    S.p=sparse(complex([0 1; 0 0]));
    S.m=sparse(complex([0 0; 1 0]));
    S.z=sparse(complex([0.5 0; 0 -0.5]));
    
elseif mult==3
    
    % Spin-one matrices are hard-coded for speed
    S.u=sparse(complex([1 0 0; 0 1 0; 0 0 1]));
    S.p=sparse(complex([0 sqrt(2) 0; 0 0 sqrt(2); 0 0 0]));
    S.m=sparse(complex([0 0 0; sqrt(2) 0 0; 0 sqrt(2) 0]));
    S.z=sparse(complex([1 0 0; 0 0 0; 0 0 -1]));
        
else
    
    % Everything else gets generated
    spin=(mult-1)/2; prjs=((mult-1):-1:0)-spin;
    S.p=spdiags(sqrt(spin*(spin+1)-prjs.*(prjs+1))',+1,mult,mult);
    S.m=spdiags(sqrt(spin*(spin+1)-prjs.*(prjs-1))',-1,mult,mult);
    S.z=spdiags(prjs',0,mult,mult); S.u=speye(mult,mult);
    S.p=complex(S.p); S.m=complex(S.m);
    S.z=complex(S.z); S.u=complex(S.u);
    
end

% X and Y operators are combinations
S.x=(S.p+S.m)/2; S.y=(S.p-S.m)/2i;
    
end

% Consistency enforcement
function grumble(mult)
if (~isnumeric(mult))||(~isreal(mult))||...
   (~isscalar(mult))||(mod(mult,1)~=0)||(mult<1)
    error('mult must be a positive real integer.');
end
end

% Any refusal to recognize reality, for any reason whatever, has dis-
% astrous consequences. There are no evil thoughts except one -- the
% refusal to think. Don't ignore your own desires... Don't sacrifice
% them. Examine their cause.
%
% Ayn Rand, "Atlas Shrugged"

