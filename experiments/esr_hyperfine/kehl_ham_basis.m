% Cached operator basis for Kehl ENDOR Hamiltonians. Syntax:
%
%      ham=kehl_ham_basis(parameters,paramsENDOR,spin_system)
%
% Parameters:
%
%   parameters       - Kehl ENDOR context parameter structure.
%   paramsENDOR      - map containing ENDOR parameters.
%   spin_system      - Zeeman-Liouville Spinach spin system.
%
% Outputs:
%
%   ham              - structure containing cached superoperators.
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

    % Cache electron and ENDOR nuclear superoperators
    electron_idx=parameters.electron_spin_idx;
    ham.Sz=operator(spin_system,'Lz',electron_idx);
    n_endor=parameters.n_endor;
    ham.Ix=cell(1,n_endor);
    ham.Iy=cell(1,n_endor);
    ham.Iz=cell(1,n_endor);
    ham.SzIx=cell(1,n_endor);
    ham.SzIy=cell(1,n_endor);
    ham.SzIz=cell(1,n_endor);
    ham.IxIx=cell(1,n_endor);
    ham.IxIy=cell(1,n_endor);
    ham.IxIz=cell(1,n_endor);
    ham.IyIx=cell(1,n_endor);
    ham.IyIy=cell(1,n_endor);
    ham.IyIz=cell(1,n_endor);
    ham.IzIx=cell(1,n_endor);
    ham.IzIy=cell(1,n_endor);
    ham.IzIz=cell(1,n_endor);
    for n=1:n_endor
        spin_idx=parameters.endor_spins(n);
        ham.Ix{n}=operator(spin_system,'Lx',spin_idx);
        ham.Iy{n}=operator(spin_system,'Ly',spin_idx);
        ham.Iz{n}=operator(spin_system,'Lz',spin_idx);
        ham.SzIx{n}=product_comm(spin_system,{'Lz','Lx'},...
            [electron_idx spin_idx]);
        ham.SzIy{n}=product_comm(spin_system,{'Lz','Ly'},...
            [electron_idx spin_idx]);
        ham.SzIz{n}=product_comm(spin_system,{'Lz','Lz'},...
            [electron_idx spin_idx]);
        ham.IxIx{n}=product_comm(spin_system,{'Lx','Lx'},...
            [spin_idx spin_idx]);
        ham.IxIy{n}=product_comm(spin_system,{'Lx','Ly'},...
            [spin_idx spin_idx]);
        ham.IxIz{n}=product_comm(spin_system,{'Lx','Lz'},...
            [spin_idx spin_idx]);
        ham.IyIx{n}=product_comm(spin_system,{'Ly','Lx'},...
            [spin_idx spin_idx]);
        ham.IyIy{n}=product_comm(spin_system,{'Ly','Ly'},...
            [spin_idx spin_idx]);
        ham.IyIz{n}=product_comm(spin_system,{'Ly','Lz'},...
            [spin_idx spin_idx]);
        ham.IzIx{n}=product_comm(spin_system,{'Lz','Lx'},...
            [spin_idx spin_idx]);
        ham.IzIy{n}=product_comm(spin_system,{'Lz','Ly'},...
            [spin_idx spin_idx]);
        ham.IzIz{n}=product_comm(spin_system,{'Lz','Lz'},...
            [spin_idx spin_idx]);
    end

    % Cache legacy first-nucleus dipolar superoperator shapes
    ham.dip_oper={};
    if parameters.dipolar_active&&(n_endor>1)
        for n=2:n_endor
            ham.dip_oper{n-1}=legacy_dip_oper(spin_system,parameters,n);
        end
    end

    % Store the basis in the process-local cache
    basis_cache(cache_key)=ham;

end

% Commutation superoperator for an ordered product of spin operators
function A=product_comm(spin_system,operators,spins)

    % Initialise left and right product superoperators
    dim=prod(spin_system.comp.mults)^2;
    left_op=speye(dim);
    right_op=speye(dim);

    % Build the left product in physical operator order
    for n=1:numel(operators)
        left_op=left_op*operator(spin_system,operators{n},spins(n),'left');
    end

    % Build the right product in reverse order
    for n=numel(operators):-1:1
        right_op=right_op*operator(spin_system,operators{n},spins(n),'right');
    end

    % Return the commutation superoperator
    A=left_op-right_op;

end

% Legacy secular dipolar superoperator shape
function D_oper=legacy_dip_oper(spin_system,parameters,spin_idx)

    % Get first and current ENDOR spin numbers
    spin_a=parameters.endor_spins(1);
    spin_b=parameters.endor_spins(spin_idx);

    % Assemble the scalar nuclear spin product superoperator
    H=product_comm(spin_system,{'Lx','Lx'},[spin_a spin_b])+...
        product_comm(spin_system,{'Lx','Ly'},[spin_a spin_b])+...
        product_comm(spin_system,{'Lx','Lz'},[spin_a spin_b])+...
        product_comm(spin_system,{'Ly','Lx'},[spin_a spin_b])+...
        product_comm(spin_system,{'Ly','Ly'},[spin_a spin_b])+...
        product_comm(spin_system,{'Ly','Lz'},[spin_a spin_b])+...
        product_comm(spin_system,{'Lz','Lx'},[spin_a spin_b])+...
        product_comm(spin_system,{'Lz','Ly'},[spin_a spin_b])+...
        product_comm(spin_system,{'Lz','Lz'},[spin_a spin_b]);

    % Return the legacy axial operator combination
    D_oper=(3*product_comm(spin_system,{'Lz','Lz'},[spin_a spin_b])-H)/2;

end

% Build a conservative cache key
function cache_key=ham_key(spin_system,parameters)

    % Identify the Spinach basis
    if isfield(spin_system.bas,'basis_hash')
        basis_id=spin_system.bas.basis_hash;
    elseif isfield(spin_system.bas,'basis')
        basis_id=spin_system.bas.basis;
    else
        basis_id={spin_system.bas.formalism,spin_system.bas.approximation};
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
    if ~strcmp(spin_system.bas.formalism,'zeeman-liouv')
        error('spin_system must use zeeman-liouv formalism.');
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

