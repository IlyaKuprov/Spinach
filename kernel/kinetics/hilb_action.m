% Applies a direct-sum Liouville map to a block-diagonal Hilbert state.
% Syntax:
%
%                  deriv=hilb_action(spin_system,L,t,rho)
%
% Parameters:
%
%    spin_system - zeeman-hilb system with compiled substance offsets
%
%    L           - sum(D_n^2)-square derivative map, or @(t,eta) returning
%                  that map; eta stacks vectorised local density matrices
%
%    t           - finite real scalar time in seconds
%
%    rho         - sum(D_n)-square block-diagonal Hilbert density matrix
%
% Outputs:
%
%    deriv       - block-diagonal matrix containing the mapped derivatives
%
% L is a derivative map, without the Hamiltonian propagation factor -1i.
% Only physical substance blocks are vectorised; inter-species coherences
% are rejected rather than discarded. Concentrations are never divided out.
%
% ilya.kuprov@weizmann.ac.il

function deriv=hilb_action(spin_system,L,t,rho)

% Check the caller-controlled matrix and time
grumble(spin_system,L,t,rho);

% Vectorise only the physical density-matrix blocks
blocks=cell(spin_system.bas.nsubst,1);
for n=1:spin_system.bas.nsubst
    idx=(spin_system.bas.offsets(n)+1):spin_system.bas.offsets(n+1);
    blocks{n}=rho(idx,idx);
end
eta=hilb2liouv(blocks,'statevec');

% Apply a fixed or time-dependent direct-sum derivative map
if isa(L,'function_handle'), L=L(t,eta); end
eta=L*eta; offsets=[0;cumsum(spin_system.bas.nstates.^2)];

% Restore the matrix blocks without inter-species coherences
for n=1:spin_system.bas.nsubst
    idx=(offsets(n)+1):offsets(n+1);
    blocks{n}=reshape(eta(idx),spin_system.bas.nstates(n),spin_system.bas.nstates(n));
end
deriv=blkdiag(blocks{:});

end

% Consistency enforcement
function grumble(spin_system,L,t,rho)
if ~strcmp(spin_system.bas.formalism,'zeeman-hilb')
    error('Spinach:hilb_action:formalism','matrix actions require zeeman-hilb formalism.');
end
if ~isnumeric(t)||~isreal(t)||~isscalar(t)||~isfinite(t)
    error('Spinach:hilb_action:time','t must be a finite real scalar.');
end
if ~isnumeric(rho)||~isequal(size(rho),[spin_system.bas.offsets(end) spin_system.bas.offsets(end)])||...
   any(~isfinite(rho),'all')
    error('Spinach:hilb_action:state','rho must be a finite matrix of the compiled Hilbert dimensions.');
end
for n=1:spin_system.bas.nsubst
    idx=(spin_system.bas.offsets(n)+1):spin_system.bas.offsets(n+1);
    if nnz(rho(idx,:))~=nnz(rho(idx,idx))
        error('Spinach:hilb_action:crossSubstance','rho must not contain inter-substance coherences.');
    end
end
dim=sum(spin_system.bas.nstates.^2);
if ~isa(L,'function_handle')&&(~isnumeric(L)||~isequal(size(L),[dim dim])||any(~isfinite(L),'all'))
    error('Spinach:hilb_action:map','L must be a finite direct-sum Liouville matrix or a function handle.');
end
end


