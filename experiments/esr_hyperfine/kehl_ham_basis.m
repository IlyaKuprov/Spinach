% Cached operator basis for Kehl ENDOR Hamiltonians. Syntax:
%
%      ham=kehl_ham_basis(parameters,paramsENDOR,spin_system)
%
% Parameters:
%
%   parameters       - Kehl ENDOR context parameter structure.
%   paramsENDOR      - map containing ENDOR parameters.
%   spin_system      - full Spinach spin system.
%
% Outputs:
%
%   ham              - structure containing cached operators.
%
% February 2024 A. Kehl (akehl@gwdg.de)
% May 2026 Spinach integration
%
% <https://spindynamics.org/wiki/index.php?title=kehl_ham_basis.m>

function ham=kehl_ham_basis(parameters,paramsENDOR,spin_system)

    % Check consistency
    grumble(parameters,paramsENDOR,spin_system);

    % Initialise the process-local Hamiltonian-basis cache
    persistent basis_cache
    if isempty(basis_cache)
        basis_cache=containers.Map;
    end

    % Build a conservative cache key
    cache_key=ham_key(spin_system,parameters);

    % Return the cached basis if it already exists
    if isKey(basis_cache,cache_key)
        ham=basis_cache(cache_key);
        return
    end

    % Cache electron and identity operators
    ham.Sz=full(operator(spin_system,'Lz',parameters.electron_spin_idx));
    ham.id=eye(size(ham.Sz));

    % Cache ENDOR nuclear operators
    n_endor=parameters.n_endor;
    ham.Ix=cell(1,n_endor);
    ham.Iy=cell(1,n_endor);
    ham.Iz=cell(1,n_endor);
    ham.spin_q=zeros(1,n_endor);
    for n=1:n_endor
        spin_idx=parameters.endor_spins(n);
        ham.Ix{n}=full(operator(spin_system,'Lx',spin_idx));
        ham.Iy{n}=full(operator(spin_system,'Ly',spin_idx));
        ham.Iz{n}=full(operator(spin_system,'Lz',spin_idx));
        ham.spin_q(n)=(spin_system.comp.mults(spin_idx)-1)/2;
    end

    % Cache legacy first-nucleus dipolar operator shapes
    ham.dip_oper={};
    if parameters.dipolar_active&&(n_endor>1)
        for n=2:n_endor
            ham.dip_oper{n-1}=legacy_dip_oper(ham,n);
        end
    end

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
function cache_key=ham_key(spin_system,parameters)

    % Identify the Spinach basis
    if isfield(spin_system.bas,'basis_hash')
        basis_id=spin_system.bas.basis_hash;
    else
        basis_id=spin_system.bas.basis;
    end

    % Hash all data that can affect the operator basis
    cache_key=md5_hash({spin_system.comp.isotopes,basis_id,...
        parameters.electron_spin_idx,parameters.endor_spins,...
        parameters.dipolar_active});

end

% Consistency enforcement
function grumble(parameters,paramsENDOR,spin_system)
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
    if (~isfield(parameters,'electron_spin_idx'))||...
            (~isnumeric(parameters.electron_spin_idx))||...
            (~isscalar(parameters.electron_spin_idx))
        error('parameters.electron_spin_idx must be a numeric scalar.');
    end
    if (~isfield(parameters,'endor_spins'))||...
            (~isnumeric(parameters.endor_spins))||...
            (~isvector(parameters.endor_spins))
        error('parameters.endor_spins must be a numeric vector.');
    end
    if (~isfield(parameters,'n_endor'))||...
            (~isnumeric(parameters.n_endor))||(~isscalar(parameters.n_endor))
        error('parameters.n_endor must be a numeric scalar.');
    end
end

