% Static ENDOR Liouvillian from direct Spinach operators. Syntax:
%
%      H=kehl_free_ham(parameters,paramsENDOR,spin_system,v_off_S,euler_angles)
%
% Parameters:
%
%   parameters       - Kehl ENDOR context parameter structure.
%   paramsENDOR      - map containing ENDOR parameters.
%   spin_system      - Zeeman-Liouville Spinach spin system.
%   v_off_S          - electron offset frequency.
%   euler_angles     - Kehl orientation angles stored by the context.
%
% Outputs:
%
%   H                - static Liouville-space commutation superoperator.
%
% February 2024 A. Kehl (akehl@gwdg.de)
% May 2026 Spinach integration
%
% <https://spindynamics.org/wiki/index.php?title=kehl_free_ham.m>

function H=kehl_free_ham(parameters,paramsENDOR,spin_system,v_off_S,euler_angles)

    % Default orientation is the input tensor frame
    if nargin<5
        euler_angles=[0 0 0];
    end

    % Check consistency
    grumble(parameters,paramsENDOR,spin_system,v_off_S,euler_angles);

    % Project tensor data onto the selected Kehl orientation
    terms=kehl_orient_terms(parameters,euler_angles);

    % Get the ENDOR Larmor frequencies and electron index
    v_L=paramsENDOR('v_L');
    electron_idx=parameters.electron_spin_idx;

    % Start with the electron frequency offset
    H=2*pi*v_off_S*operator(spin_system,'Lz',electron_idx);

    % Add all ENDOR nuclear terms in the legacy Kehl form
    for n=1:parameters.n_endor
        spin_idx=parameters.endor_spins(n);
        Iz=operator(spin_system,'Lz',spin_idx);
        H=H-2*pi*v_L(n)*Iz+2*pi*v_L(n)*terms.cs_zz(n)*Iz;
        H=H+2*pi*terms.hf_zz(n)*product_comm(spin_system,...
            {'Lz','Lz'},[electron_idx spin_idx])+...
            2*pi*terms.hf_zy(n)*product_comm(spin_system,...
            {'Lz','Ly'},[electron_idx spin_idx])+...
            2*pi*terms.hf_zx(n)*product_comm(spin_system,...
            {'Lz','Lx'},[electron_idx spin_idx]);
        if parameters.nqi_active
            if parameters.Bterm
                H=H+nqi_full(spin_system,terms.nqi(n,:,:),spin_idx);
            else
                H=H+3*pi*terms.nqi_zz(n)*product_comm(spin_system,...
                    {'Lz','Lz'},[spin_idx spin_idx]);
            end
        end
    end

    % Add the legacy first-nucleus dipolar correction
    if parameters.dipolar_active
        for n=1:numel(terms.dip_zz)
            H=H+2*pi*terms.dip_zz(n)*legacy_dip_oper(...
                spin_system,parameters,n+1);
        end
    end

end

% Full quadrupolar B-term contribution
function H=nqi_full(spin_system,nqi_tensor,spin_idx)

    % Reshape the selected tensor slice
    Q=squeeze(nqi_tensor);

    % Assemble the full quadrupolar Liouvillian
    H=Q(1,1)*product_comm(spin_system,{'Lx','Lx'},[spin_idx spin_idx])+...
        Q(1,2)*product_comm(spin_system,{'Lx','Ly'},[spin_idx spin_idx])+...
        Q(1,3)*product_comm(spin_system,{'Lx','Lz'},[spin_idx spin_idx])+...
        Q(2,1)*product_comm(spin_system,{'Ly','Lx'},[spin_idx spin_idx])+...
        Q(2,2)*product_comm(spin_system,{'Ly','Ly'},[spin_idx spin_idx])+...
        Q(2,3)*product_comm(spin_system,{'Ly','Lz'},[spin_idx spin_idx])+...
        Q(3,1)*product_comm(spin_system,{'Lz','Lx'},[spin_idx spin_idx])+...
        Q(3,2)*product_comm(spin_system,{'Lz','Ly'},[spin_idx spin_idx])+...
        Q(3,3)*product_comm(spin_system,{'Lz','Lz'},[spin_idx spin_idx]);

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

% Consistency enforcement
function grumble(parameters,paramsENDOR,spin_system,v_off_S,euler_angles)
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
    if ~isnumeric(v_off_S)
        error('v_off_S must be numeric.');
    end
    if (~isnumeric(euler_angles))||(~isvector(euler_angles))||...
            (numel(euler_angles)~=3)
        error('euler_angles must be a three-element numeric vector.');
    end
end

