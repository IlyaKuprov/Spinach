% Singlet-singlet RYDMR experiment using the full kinetics superoper-
% ator - computes the singlet yield of a radical pair recombination 
% reaction. Syntax:
%
%             A=rydmr(spin_system,parameters,H,R,K)
%
% where H is the Hamiltonian commutation superoperator in zero ex-
% ternal field, R is the relaxation superoperator and K is the che-
% mical kinetics superoperator. Parameters:
%
%     parameters.tol     -  BICG solver tolerance,
%                           1e-2 is generally good
%
% Outputs:
%
%     A - fractional singlet yield
%
% ilya.kuprov@weizmann.ac.uk
% h.j.hogben@chem.ox.ac.uk
% peter.hore@chem.ox.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=rydmr.m>

function A=rydmr(spin_system,parameters,H,R,K)

% Check consistency
grumble(spin_system,parameters,H,R,K);

% Get the two-electron singlet state
S=singlet(spin_system,spin_system.chem.rp_electrons(1),...
                      spin_system.chem.rp_electrons(2));

% Compose Liouvillian
L=H+1i*R+1i*K;
                  
% Normalize the singlet
S=S/norm(S,2);

% Move to GPU if needed
if ismember('gpu',spin_system.sys.enable)
    L=gpuArray(L); S=gpuArray(S);
end

% Compute singlet yield
A=spin_system.chem.rp_rates(1)*...
  imag(S'*bicg(L,S,parameters.tol,numel(S)));

% Gather from GPU if needed
if ismember('gpu',spin_system.sys.enable)
    A=gather(A);
end
    
end

% Consistency enforcement
function grumble(spin_system,parameters,H,R,K) %#ok<INUSL>
if (~isnumeric(H))||(~isnumeric(R))||(~isnumeric(K))||...
   (~ismatrix(H))||(~ismatrix(R))||(~ismatrix(K))
    error('H, R and K arguments must be matrices.');
end
if (~all(size(H)==size(R)))||(~all(size(R)==size(K)))
    error('H, R and K matrices must have the same dimension.');
end
end

% You can stand on the shoulders of giants, or a big 
% enough pile of dwarfs, works either way.
%
% Internet folklore

