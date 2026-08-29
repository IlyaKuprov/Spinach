% Coherence order selection function - keeps only the specified orders
% of coherence in the state vector. This is useful as an analytical re-
% placement for complicated phase cycles. Syntax:
%
%                   rho=coherence(spin_system,rho,spec)
%
% Arguments:
%
%     rho    -  a state vector or a horizontal stack thereof;
%               in zeeman-hilb, a density matrix or a horizon-
%               tal stack thereof
%
%     spec   -  a cell array containing the specification of
%               which coherences to keep on which spins. For
%               example
%                         {{'13C',[1 -1]},{'1H',-1}} 
%
%               keeps the states that have coherence order 
%
%                      ((1 OR -1 on 13C) AND (-1 on 1H))
%
%               instead of specific spins, it is possible to
%               specify 'electrons', 'nuclei', and 'all'
%
% Outputs:
%
%   rho     - the state vector with the undesired orders of
%             spin correlations zeroed out
%
% Note: this function requires sphten-liouv, zeeman-liouv, or zeeman-
%       hilb formalism; Fokker-Planck direct products are supported
%       in the Liouville space formalisms. In zeeman-hilb, the densi-
%       ty matrices are stretched into Liouville space, filtered the-
%       re, and folded back.
%
% ilya.kuprov@weizmann.ac.il
% ledwards@cbs.mpg.de
%
% <https://spindynamics.org/wiki/index.php?title=coherence.m>

function rho=coherence(spin_system,rho,spec)

% Check consistency
grumble(spin_system,rho,spec);

% Store dimension statistics
spn_dim=size(spin_system.bas.basis,1);
if strcmp(spin_system.bas.formalism,'zeeman-hilb')
    spn_dim=spn_dim^2;
end
spc_dim=numel(rho)/spn_dim;
problem_dims=size(rho);

% Fold indirect dimensions
rho=reshape(rho,[spn_dim spc_dim]);

% Compute coherence order bookkeeping array
switch spin_system.bas.formalism

    case 'sphten-liouv'

        % Projection quantum numbers of basis states
        [~,M]=lin2lm(spin_system.bas.basis);

    case 'zeeman-liouv'

        % Projection quantum numbers of ket and bra indices
        nspins=spin_system.comp.nspins;
        spns=(spin_system.comp.mults-1)/2;
        M_ket=spns-spin_system.bas.basis(:,1:nspins)+1;
        M_bra=spns-spin_system.bas.basis(:,(nspins+1):end)+1;

        % Coherence orders of stretched density matrix elements
        M=M_ket-M_bra;

    case 'zeeman-hilb'

        % Projection quantum numbers of the Zeeman basis
        hdim=size(spin_system.bas.basis,1);
        spns=(spin_system.comp.mults-1)/2;
        M_lvl=spns-spin_system.bas.basis+1;

        % Coherence orders of stretched density matrix elements
        M=repmat(M_lvl,[hdim 1])-kron(M_lvl,ones(hdim,1));

end

% Preallocate state mask array
state_mask=false(spn_dim,numel(spec));

% Loop over specifications
for n=1:numel(spec)

    % Parse spin specification
    if ischar(spec{n}{1})
        
        % Symbolic specification
        if strcmp(spec{n}{1},'all')
            spins=1:numel(spin_system.comp.isotopes);
        elseif strcmp(spec{n}{1},'electrons')
            spins=cellfun(@iselectron,spin_system.comp.isotopes);
        elseif strcmp(spec{n}{1},'nuclei')
            spins=cellfun(@isnucleus,spin_system.comp.isotopes);
        else
            spins=strcmp(spec{n}{1},spin_system.comp.isotopes);
        end
        
    else
        
        % Specification by number
        spins=spec{n}{1};
        
    end
    
    % Determine coherence order of each basis state
    coherence_orders_present=sum(M(:,spins),2);
  
    % Wipe all coherence orders except those specified by the user
    state_mask(:,n)=ismember(coherence_orders_present,spec{n}{2});
        
end

% Intersect state masks
state_mask=all(state_mask,2);

% Apply the state mask
rho(~state_mask,:)=0;

% Unfold indirect dimensions
rho=reshape(rho,problem_dims);

% Report overly destructive calls
if norm(rho,1)<1e-10
    report(spin_system,'WARNING - all magnetization appears to have been destroyed by this call.');
end

end

% Consistency enforcement
function grumble(spin_system,rho,spec)
if ~ismember(spin_system.bas.formalism,{'sphten-liouv','zeeman-liouv','zeeman-hilb'})
    error('analytical coherence order selection is only available for sphten-liouv, zeeman-liouv, and zeeman-hilb formalisms.');
end
if ~isnumeric(rho)
    error('the state vector(s) must be numeric.');
end
if mod(numel(rho),size(spin_system.bas.basis,1))~=0
    error('the number of elements in rho must be a multiple of the dimension of the spin state space.');
end
if ~iscell(spec)
    error('spec must be a cell array of cell arrays, e.g. {{''13C'',[1 -1]},{''1H'',-1}}');
end
for n=1:numel(spec)
    if (~iscell(spec{n}))||(numel(spec{n})~=2)
        error('each element of spec must be a two-element cell array.');
    end
    if ischar(spec{n}{1})
        if (~ismember(spec{n}{1},{'all','electrons','nuclei'}))&&...
           (~ismember(spec{n}{1},spin_system.comp.isotopes))
            error('spec spin string must be ''all'', ''electrons'', ''nuclei'', or an isotope present in the system.');
        end
    elseif isnumeric(spec{n}{1})
        if (~isvector(spec{n}{1}))||(~isreal(spec{n}{1}))||any(spec{n}{1}<1)||...
           any(mod(spec{n}{1},1)~=0)||any(spec{n}{1}>spin_system.comp.nspins)
            error('spin numbers in spec must be a vector of positive integers within the system bounds.');
        end
    else
        error('spin specification in spec must be ''all'', ''electrons'', ''nuclei'', an isotope string, or a vector of spin numbers.');
    end
    if (~isnumeric(spec{n}{2}))||(~isvector(spec{n}{2}))||(~isreal(spec{n}{2}))||...
       any(mod(spec{n}{2},1)~=0)
        error('coherence orders in spec must be a vector of real integers.');
    end
end
end

% There was a time when men were afraid that somebody would reveal
% some secret of theirs that was unknown to their fellows. Nowadays,
% they're afraid that somebody will name what everybody knows. Have
% you practical people ever thought that that's all it would take
% to blast your whole, big, complex structure, with all your laws
% and guns -- just somebody naming the exact nature of what it is
% you're doing?
%
% Ayn Rand, "Atlas Shrugged"

