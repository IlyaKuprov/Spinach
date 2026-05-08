% Cached operators for Kehl ENDOR kernels. Syntax:
%
%      ops=kehl_operator_basis(spin_system,n_endor,n_spin_systems)
%
% Parameters:
%
%   spin_system      - reduced Spinach spin system.
%   n_endor          - number of ENDOR nuclei.
%   n_spin_systems   - number of separated spin systems.
%
% Outputs:
%
%   ops              - structure containing cached operators and states.
%
% February 2024 A. Kehl (akehl@gwdg.de)
% May 2026 Spinach integration
%
% <https://spindynamics.org/wiki/index.php?title=kehl_operator_basis.m>

function ops=kehl_operator_basis(spin_system,n_endor,n_spin_systems)

    % Default is the compact ENDOR spin system
    if nargin<3
        n_spin_systems=1;
    end

    % Check consistency
    grumble(spin_system,n_endor,n_spin_systems);

    % Initialise the process-local operator cache
    persistent operator_cache
    if isempty(operator_cache)
        operator_cache=containers.Map;
    end

    % Build a conservative cache key
    cache_key=operator_key(spin_system,n_endor,n_spin_systems);

    % Return cached operators when available
    if isKey(operator_cache,cache_key)
        ops=operator_cache(cache_key);
        return
    end

    % Build electron operators once
    ops.Sx=full(operator(spin_system,'Lx',1));
    ops.Sy=full(operator(spin_system,'Ly',1));
    ops.Sz=full(operator(spin_system,'Lz',1));

    % Build nuclear operators once
    n_nuc=spin_system.comp.nspins-1;
    ops.Ix=cell(1,n_nuc);
    ops.Iy=cell(1,n_nuc);
    ops.Iz=cell(1,n_nuc);
    for n=1:n_nuc
        spin_idx=n+1;
        ops.Ix{n}=full(operator(spin_system,'Lx',spin_idx));
        ops.Iy{n}=full(operator(spin_system,'Ly',spin_idx));
        ops.Iz{n}=full(operator(spin_system,'Lz',spin_idx));
    end

    % Reuse one nuclear operator set in separated-spin-system mode
    if n_spin_systems>1
        for n=2:n_endor
            ops.Ix{n}=ops.Ix{1};
            ops.Iy{n}=ops.Iy{1};
            ops.Iz{n}=ops.Iz{1};
        end
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
    if n_spin_systems>1
        ops.Ix_rf=ops.Ix{1};
        ops.Iy_rf=ops.Iy{1};
        ops.Iz_rf=ops.Iz{1};
    else
        for n=1:n_endor
            ops.Ix_rf=ops.Ix_rf+ops.Ix{n};
            ops.Iy_rf=ops.Iy_rf+ops.Iy{n};
            ops.Iz_rf=ops.Iz_rf+ops.Iz{n};
        end
    end

    % Store dimensions and identity for small repeated constructions
    ops.dim=size(ops.Sz,1);
    ops.eye=eye(ops.dim);

    % Store the operators in the process-local cache
    operator_cache(cache_key)=ops;

end

% Build a conservative operator cache key
function cache_key=operator_key(spin_system,n_endor,n_spin_systems)

    % Identify the reduced Spinach basis
    if isfield(spin_system.bas,'basis_hash')
        basis_id=spin_system.bas.basis_hash;
    else
        basis_id=spin_system.bas.basis;
    end

    % Hash all data that can affect the operator basis
    cache_key=md5_hash({spin_system.comp.isotopes,basis_id,n_endor,...
        n_spin_systems});

end

% Consistency enforcement
function grumble(spin_system,n_endor,n_spin_systems)
    if (~isstruct(spin_system))||(~isfield(spin_system,'bas'))||...
            (~isfield(spin_system,'comp'))
        error('spin_system must be a Spinach spin system structure.');
    end
    if (~isnumeric(n_endor))||(~isscalar(n_endor))||...
            (n_endor<0)||mod(n_endor,1)~=0
        error('n_endor must be a non-negative integer.');
    end
    if (~isnumeric(n_spin_systems))||(~isscalar(n_spin_systems))||...
            (n_spin_systems<1)||mod(n_spin_systems,1)~=0
        error('n_spin_systems must be a positive integer.');
    end
    if (n_spin_systems>1)&&(spin_system.comp.nspins<2)
        error('separated spin systems require at least one nuclear spin.');
    end
end

