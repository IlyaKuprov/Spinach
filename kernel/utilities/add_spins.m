% Reduction of direct products of two su(2) irreps. Syntax:
%
%             [mult,proj]=add_spins(spin_a,spin_b)
%
% Parameters:
%
%     spin_a  - quantum number of the first spin,
%               an integer or a half-integer
%
%     spin_b  - quantum number of the second spin,
%               an integer or a half-integer
%
% Outputs:
%
%     mult    - multiplicities corresponding to
%               the values of the total spin that
%               are present
%
%     proj    - projectors that reduce the direct
%               product representation, a cell ar-
%               ray of matrices
%
% ilya.kuprov@weizmann.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=add_spins.m>

function [mult,proj]=add_spins(spin_a,spin_b)

% Check consistency
grumble(spin_a,spin_b);

% Get the individual irreps
spin_a=pauli(2*spin_a+1); 
spin_b=pauli(2*spin_b+1);

% Build direct product rep generators
Sx=kron(spin_a.u,spin_b.x)+kron(spin_a.x,spin_b.u);
Sy=kron(spin_a.u,spin_b.y)+kron(spin_a.y,spin_b.u);
Sz=kron(spin_a.u,spin_b.z)+kron(spin_a.z,spin_b.u);

% Diagonalise Casimir operator
[V,D]=eig(full(Sx^2+Sy^2+Sz^2),'vector');

% Index Casimir operator eigenvalues
[S_sq,~,idx_b]=unique(uint32(D));

% Make projectors
proj=cell(1,numel(S_sq));
for n=1:numel(S_sq)
    proj{n}=V(:,idx_b==n);
end

% Canonicalise projectors
for n=1:numel(S_sq)
    
    % Set the expectations
    mult(n)=size(proj{n},2); %#ok<AGROW>
    Q=pauli(mult(n));
    
    % Set Sz eigenvalue order
    Sz_irr=proj{n}'*Sz*proj{n};
    [V,D]=eig(full(Sz_irr),'vector');
    [~,idx]=sort(D,'descend'); 
    
    % Require diagonal Sz
    V=V(:,idx); proj{n}=proj{n}*V;
    if norm(proj{n}'*Sz*proj{n}-Q.z,1)>sqrt(eps)
        error('irrep canonicalisation failed.');
    end
    
    % Require real and positive Sx
    while any(proj{n}'*Sx*proj{n}+sqrt(eps)<0,'all')
        Sx_irr=proj{n}'*Sx*proj{n};
        for k=1:size(Sx_irr,1)
            if any(Sx_irr(:,k)+sqrt(eps)<0)
                proj{n}(:,k)=-proj{n}(:,k); break;
            end
        end
    end
    if (norm(proj{n}'*Sx*proj{n}-Q.x,1)>sqrt(eps))||...
       (norm(proj{n}'*Sy*proj{n}-Q.y,1)>sqrt(eps))
        error('irrep canonicalisation failed.');
    end
    
end

end

% Consistency enforcement
function grumble(spin_a,spin_b)
if (~isnumeric(spin_a))||(~isreal(spin_a))||...
   (~isscalar(spin_a))||(mod(2*spin_a+1,1)~=0)||...
   (spin_a<1/2)
    error('spin_a must be a positive integer or half-integer.');
end
if (~isnumeric(spin_b))||(~isreal(spin_b))||...
   (~isscalar(spin_b))||(mod(2*spin_b+1,1)~=0)||...
   (spin_b<1/2)
    error('spin_b must be a positive integer or half-integer.');
end
end

% You must gather your party before venturing forth.
%
% Narrator in Baldur's Gate

