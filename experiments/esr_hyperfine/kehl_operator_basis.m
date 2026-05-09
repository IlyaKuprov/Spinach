% Cached operators for Kehl ENDOR kernels. Syntax:
%
%      ops=kehl_operator_basis(spin_system,parameters)
%
% Parameters:
%
%   spin_system      - full Spinach spin system.
%   parameters       - Kehl ENDOR context parameter structure.
%
% Outputs:
%
%   ops              - structure containing cached operators and states.
%
% February 2024 A. Kehl (akehl@gwdg.de)
% May 2026 Spinach integration
%
% <https://spindynamics.org/wiki/index.php?title=kehl_operator_basis.m>

function ops=kehl_operator_basis(spin_system,parameters)

    % Check consistency
    grumble(spin_system,parameters);

    % Initialise the process-local operator cache
    persistent operator_cache
    if isempty(operator_cache)
        operator_cache=containers.Map;
    end

    % Build a conservative cache key
    cache_key=operator_key(spin_system,parameters);

    % Return cached operators when available
    if isKey(operator_cache,cache_key)
        ops=operator_cache(cache_key);
        return
    end

    % Build electron operators once
    electron_idx=parameters.electron_spin_idx;
    ops.Sx=full(operator(spin_system,'Lx',electron_idx));
    ops.Sy=full(operator(spin_system,'Ly',electron_idx));
    ops.Sz=full(operator(spin_system,'Lz',electron_idx));

    % Build ENDOR nuclear operators once
    n_endor=parameters.n_endor;
    ops.Ix=cell(1,n_endor);
    ops.Iy=cell(1,n_endor);
    ops.Iz=cell(1,n_endor);
    for n=1:n_endor
        spin_idx=parameters.endor_spins(n);
        ops.Ix{n}=full(operator(spin_system,'Lx',spin_idx));
        ops.Iy{n}=full(operator(spin_system,'Ly',spin_idx));
        ops.Iz{n}=full(operator(spin_system,'Lz',spin_idx));
    end

    % Precompute legacy diagonal-frame relaxation operators
    [ops.Sx_D,ops.Sy_D,ops.Sz_D,ops.Ix_D,ops.Iy_D,ops.Iz_D]=...
        kehl_diag_ops(ops.Sx,ops.Sy,ops.Sz,ops.Ix,ops.Iy,ops.Iz,n_endor);

    % Precompute the legacy density matrix without normalisation
    ops.rho_z=-ops.Sz;

    % Precompute the repeated RF correction operators
    ops.Ix_rf=zeros(size(ops.Sz));
    ops.Iy_rf=zeros(size(ops.Sz));
    ops.Iz_rf=zeros(size(ops.Sz));
    for n=1:n_endor
        ops.Ix_rf=ops.Ix_rf+ops.Ix{n};
        ops.Iy_rf=ops.Iy_rf+ops.Iy{n};
        ops.Iz_rf=ops.Iz_rf+ops.Iz{n};
    end

    % Store dimensions and identity for small repeated constructions
    ops.dim=size(ops.Sz,1);
    ops.eye=eye(ops.dim);

    % Store the operators in the process-local cache
    operator_cache(cache_key)=ops;

end

% Build a conservative operator cache key
function cache_key=operator_key(spin_system,parameters)

    % Identify the Spinach basis
    if isfield(spin_system.bas,'basis_hash')
        basis_id=spin_system.bas.basis_hash;
    else
        basis_id=spin_system.bas.basis;
    end

    % Hash all data that can affect the operator basis
    cache_key=md5_hash({spin_system.comp.isotopes,basis_id,...
        parameters.electron_spin_idx,parameters.endor_spins});

end

% Consistency enforcement
function grumble(spin_system,parameters)
    if (~isstruct(spin_system))||(~isfield(spin_system,'bas'))||...
            (~isfield(spin_system,'comp'))
        error('spin_system must be a Spinach spin system structure.');
    end
    if ~isstruct(parameters)
        error('parameters must be a structure.');
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

