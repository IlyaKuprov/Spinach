% Builds a transition support projector for STMAS in a complete
% Zeeman Liouville basis. Syntax:
%
%             P=stmas_filter(spin_system,spin_idx,level_pairs)
%
% Parameters:
%
%     spin_system - single-substance system with a zeeman-liouv basis
%     spin_idx    - index of the selected spin in sys.isotopes
%     level_pairs - N-by-2 array of ordered (row, column) Zeeman level
%                   indices, numbered from highest to lowest m
%
% Outputs:
%
%     P           - sparse projector on column-vectorised spin states
%
% Each listed pair retains its whole spectator block without weighting
% by a transition moment. Repeated pairs count once; a 0-by-2 list
% retains nothing. The reverse coherence requires its own ordered pair.
%
% talos@spindynamics.org

function P=stmas_filter(spin_system,spin_idx,level_pairs)

% Check consistency
grumble(spin_system,spin_idx,level_pairs);

% Map product-basis indices to levels of the selected spin
spin_idx=double(spin_idx);
mults=spin_system.comp.mults;
mult=mults(spin_idx);
levels=kron(ones(prod(mults(1:spin_idx-1)),1),...
            kron((1:mult)',ones(prod(mults(spin_idx+1:end)),1)));

% Mark selected density matrix elements and project their linear indices
level_pairs=double(level_pairs);
support=sparse(level_pairs(:,1),level_pairs(:,2),true(size(level_pairs,1),1),mult,mult);
keep=find(support(levels,levels));
liouv_dim=spin_system.bas.offsets(end);
P=sparse(keep,keep,ones(size(keep)),liouv_dim,liouv_dim);

end

% Consistency enforcement
function grumble(spin_system,spin_idx,level_pairs)
if (~isstruct(spin_system))||(~isfield(spin_system,'bas'))||...
   (~strcmp(spin_system.bas.formalism,'zeeman-liouv'))||(spin_system.bas.nsubst~=1)
    error('Spinach:stmas_filter:basis','a single-substance zeeman-liouv basis is required.');
end
if (~isnumeric(spin_idx))||(~isscalar(spin_idx))||(~isreal(spin_idx))||...
   (mod(spin_idx,1)~=0)||(spin_idx<1)||(spin_idx>spin_system.comp.nspins)||...
   (~strcmp(spin_system.comp.types{spin_idx},'S'))
    error('Spinach:stmas_filter:spin','spin_idx must index a spin in the system.');
end
if (~isnumeric(level_pairs))||(~ismatrix(level_pairs))||(size(level_pairs,2)~=2)||...
   (~isreal(level_pairs))||any(mod(level_pairs,1)~=0,'all')||...
   any(level_pairs<1,'all')||any(level_pairs>spin_system.comp.mults(spin_idx),'all')
    error('Spinach:stmas_filter:levels','level_pairs must be an N-by-2 array of valid Zeeman level indices.');
end
end


