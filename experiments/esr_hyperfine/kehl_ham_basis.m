% Cached operator basis for Kehl ENDOR Hamiltonians. Syntax:
%
%      ham=kehl_ham_basis(parameters,paramsENDOR,spin_system,...
%                         spin_map,use_dipolar,term_map)
%
% Parameters:
%
%   parameters       - Kehl ENDOR context parameter structure.
%   paramsENDOR      - map containing ENDOR parameters.
%   spin_system      - reduced Spinach spin system.
%   spin_map         - ENDOR data indices represented by the reduced system.
%   use_dipolar      - flag enabling the legacy first-nucleus dipolar term.
%   term_map         - optional interaction-index override map.
%
% Outputs:
%
%   ham              - structure containing cached operators and term maps.
%
% February 2024 A. Kehl (akehl@gwdg.de)
% May 2026 Spinach integration
%
% <https://spindynamics.org/wiki/index.php?title=kehl_ham_basis.m>

function ham=kehl_ham_basis(parameters,paramsENDOR,spin_system,...
        spin_map,use_dipolar,term_map)

    % Default is to preserve the legacy dipolar correction
    if nargin<5
        use_dipolar=true;
    end

    % Default is to use the same data index for every interaction type
    if nargin<6
        term_map=struct();
    end

    % Check consistency
    grumble(parameters,paramsENDOR,spin_system,spin_map,use_dipolar,term_map);

    % Initialise the process-local Hamiltonian-basis cache
    persistent basis_cache
    if isempty(basis_cache)
        basis_cache=containers.Map;
    end

    % Normalise the spin-map and term-map shapes
    spin_map=spin_map(:).';
    term_map=complete_map(term_map,spin_map);

    % Build a conservative cache key
    cache_key=ham_key(spin_system,spin_map,use_dipolar,term_map);

    % Return the cached basis if it already exists
    if isKey(basis_cache,cache_key)
        ham=basis_cache(cache_key);
        return
    end

    % Cache electron and identity operators
    ham.Sz=full(operator(spin_system,'Lz',1));
    ham.id=eye(size(ham.Sz));

    % Cache local nuclear operators
    n_nuc=numel(spin_map);
    for n=1:n_nuc
        spin_idx=n+1;
        ham.Ix{n}=full(operator(spin_system,'Lx',spin_idx));
        ham.Iy{n}=full(operator(spin_system,'Ly',spin_idx));
        ham.Iz{n}=full(operator(spin_system,'Lz',spin_idx));
        ham.spin_q(n)=(spin_system.comp.mults(spin_idx)-1)/2;
    end

    % Cache legacy first-nucleus dipolar operator shapes
    ham.dip_oper={};
    if use_dipolar&&parameters.dipolar_active&&(n_nuc>1)
        for n=2:n_nuc
            ham.dip_oper{n-1}=legacy_dip_oper(ham,n);
        end
    end

    % Store index maps needed by the Hamiltonian assembler
    ham.spin_map=spin_map;
    ham.term_map=term_map;
    ham.use_dipolar=logical(use_dipolar);

    % Store the basis in the process-local cache
    basis_cache(cache_key)=ham;

end

% Legacy secular dipolar operator shape
function D_oper=legacy_dip_oper(ham,spin_idx)

    % Assemble the scalar nuclear spin product
    H=ham.Ix{1}*ham.Ix{spin_idx}+ham.Ix{1}*ham.Iy{spin_idx}+...
        ham.Ix{1}*ham.Iz{spin_idx}+ham.Iy{1}*ham.Ix{spin_idx}+...
        ham.Iy{1}*ham.Iy{spin_idx}+ham.Iy{1}*ham.Iz{spin_idx}+...
        ham.Iz{1}*ham.Ix{spin_idx}+ham.Iz{1}*ham.Iy{spin_idx}+...
        ham.Iz{1}*ham.Iz{spin_idx};

    % Return the legacy axial operator combination
    D_oper=(3*ham.Iz{1}*ham.Iz{spin_idx}-H)/2;

end

% Build a conservative cache key
function cache_key=ham_key(spin_system,spin_map,use_dipolar,term_map)

    % Identify the reduced Spinach basis
    if isfield(spin_system.bas,'basis_hash')
        basis_id=spin_system.bas.basis_hash;
    else
        basis_id=spin_system.bas.basis;
    end

    % Hash all data that can affect the operator basis
    cache_key=md5_hash({spin_system.comp.isotopes,basis_id,spin_map,...
        logical(use_dipolar),term_map});

end

% Complete omitted legacy term-map fields
function term_map=complete_map(term_map,spin_map)
    if ~isfield(term_map,'zeeman')
        term_map.zeeman=spin_map;
    end
    if ~isfield(term_map,'hyperfine')
        term_map.hyperfine=spin_map;
    end
    if ~isfield(term_map,'nqi')
        term_map.nqi=spin_map;
    end
    if ~isfield(term_map,'cs')
        term_map.cs=spin_map;
    end
    if ~isfield(term_map,'cs_larmor')
        term_map.cs_larmor=term_map.zeeman;
    end
    term_map.zeeman=term_map.zeeman(:).';
    term_map.hyperfine=term_map.hyperfine(:).';
    term_map.nqi=term_map.nqi(:).';
    term_map.cs=term_map.cs(:).';
    term_map.cs_larmor=term_map.cs_larmor(:).';
end

% Consistency enforcement
function grumble(parameters,paramsENDOR,spin_system,spin_map,...
        use_dipolar,term_map)
    if ~isstruct(parameters)
        error('parameters must be a structure.');
    end
    if ~isa(paramsENDOR,'containers.Map')
        error('paramsENDOR must be a containers.Map object.');
    end
    if (~isstruct(spin_system))||(~isfield(spin_system,'bas'))||...
            (~isfield(spin_system,'comp'))
        error('spin_system must be a Spinach spin system structure.');
    end
    if (~isnumeric(spin_map))||(~isvector(spin_map))||isempty(spin_map)||...
            any(spin_map<1)||any(mod(spin_map,1)~=0)
        error('spin_map must be a vector of positive integers.');
    end
    if spin_system.comp.nspins~=numel(spin_map)+1
        error('spin_map must match the number of nuclear spins in spin_system.');
    end
    if (~islogical(use_dipolar))&&(~isnumeric(use_dipolar))
        error('use_dipolar must be logical or numeric.');
    end
    if ~isstruct(term_map)
        error('term_map must be a structure.');
    end
end

