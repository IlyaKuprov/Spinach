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
% Chemistry must contain exactly one named singlet-selector record
% with a numeric rate; its electron indices define the initial pair.
% Both Haberkorn and Jones-Hore singlet selectors are accepted. All
% reaction records must have empty products: this full-space resolvent
% requires untracked loss, not population stored in a stationary product.
%
% Outputs:
%
%     A - fractional singlet yield
%
% ilya.kuprov@weizmann.ac.il
% h.j.hogben@chem.ox.ac.uk
% peter.hore@chem.ox.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=rydmr.m>

function A=rydmr(spin_system,parameters,H,R,K)

% Check consistency
grumble(spin_system,parameters,H,R,K);

% Locate the singlet reaction channel
channels=cellfun(@(r)isfield(r,'selector')&&ischar(r.selector{1})&&...
                ismember(r.selector{1},{'singlet','jones-hore-singlet'}),...
                spin_system.chem.reactions);
reaction=spin_system.chem.reactions{channels};

% Get the two-electron singlet state
S=singlet(spin_system,reaction.selector{2}(1),reaction.selector{2}(2));

% Compose Liouvillian
L=H+1i*R+1i*K;
                  
% Normalize the singlet
S=S/norm(S,2);

% Move to GPU if needed
if ismember('gpu',spin_system.sys.enable)
    L=gpuArray(L); S=gpuArray(S);
end

% Compute singlet yield
A=reaction.rate*...
  imag(S'*bicg(L,S,parameters.tol,numel(S)));

% Gather from GPU if needed
if ismember('gpu',spin_system.sys.enable)
    A=gather(A);
end
    
end

% Consistency enforcement
function grumble(spin_system,parameters,H,R,K)
channels=cellfun(@(r)isfield(r,'selector')&&ischar(r.selector{1})&&...
                ismember(r.selector{1},{'singlet','jones-hore-singlet'}),...
                spin_system.chem.reactions);
if nnz(channels)~=1
    error('exactly one named singlet-selector reaction is required.');
end
if any(cellfun(@(r)~isempty(r.products),spin_system.chem.reactions))
    error('Spinach:rydmr:trackedProducts',...
          'rydmr requires empty reaction products; propagate tracked products in the time domain.');
end
if ~isnumeric(spin_system.chem.reactions{channels}.rate)
    error('the singlet reaction rate must be numeric.');
end
if (~isnumeric(H))||(~isnumeric(R))||(~isnumeric(K))||...
   (~ismatrix(H))||(~ismatrix(R))||(~ismatrix(K))
    error('H, R and K arguments must be matrices.');
end
if (~all(size(H)==size(R)))||(~all(size(R)==size(K)))
    error('H, R and K matrices must have the same dimension.');
end
if ~isfield(parameters,'tol')
    error('solver tolerance should be specified in parameters.tol variable.');
end
if (~isnumeric(parameters.tol))||(~isreal(parameters.tol))||...
   (~isscalar(parameters.tol))||(parameters.tol<=0)
    error('parameters.tol must be a positive real scalar.');
end
end

% You can stand on the shoulders of giants, or a big 
% enough pile of dwarfs, works either way.
%
% Internet folklore

