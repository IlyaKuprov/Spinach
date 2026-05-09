% Static ENDOR Liouvillian from cached Spinach operators. Syntax:
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

    % Retrieve cached Cartesian superoperators
    ham=kehl_ham_basis(parameters,paramsENDOR,spin_system);

    % Project tensor data onto the selected Kehl orientation
    terms=kehl_orient_terms(parameters,euler_angles);

    % Get the ENDOR Larmor frequencies
    v_L=paramsENDOR('v_L');

    % Start with the electron frequency offset
    H=2*pi*v_off_S*ham.Sz;

    % Add all ENDOR nuclear terms in the legacy Kehl form
    for n=1:parameters.n_endor
        H=H-2*pi*v_L(n)*ham.Iz{n}+2*pi*v_L(n)*...
            terms.cs_zz(n)*ham.Iz{n};
        H=H+2*pi*terms.hf_zz(n)*ham.SzIz{n}+...
            2*pi*terms.hf_zy(n)*ham.SzIy{n}+...
            2*pi*terms.hf_zx(n)*ham.SzIx{n};
        if parameters.nqi_active
            if parameters.Bterm
                H=H+nqi_full(terms.nqi(n,:,:),ham,n);
            else
                H=H+3*pi*terms.nqi_zz(n)*ham.IzIz{n};
            end
        end
    end

    % Add the legacy first-nucleus dipolar correction
    if ~isempty(ham.dip_oper)
        for n=1:numel(ham.dip_oper)
            H=H+2*pi*terms.dip_zz(n)*ham.dip_oper{n};
        end
    end

end

% Full quadrupolar B-term contribution
function H=nqi_full(nqi_tensor,ham,spin_idx)

    % Reshape the selected tensor slice
    Q=squeeze(nqi_tensor);

    % Assemble the full quadrupolar Liouvillian
    H=Q(1,1)*ham.IxIx{spin_idx}+Q(1,2)*ham.IxIy{spin_idx}+...
        Q(1,3)*ham.IxIz{spin_idx}+Q(2,1)*ham.IyIx{spin_idx}+...
        Q(2,2)*ham.IyIy{spin_idx}+Q(2,3)*ham.IyIz{spin_idx}+...
        Q(3,1)*ham.IzIx{spin_idx}+Q(3,2)*ham.IzIy{spin_idx}+...
        Q(3,3)*ham.IzIz{spin_idx};

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

