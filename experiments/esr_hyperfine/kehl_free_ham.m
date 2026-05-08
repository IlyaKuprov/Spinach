% Static ENDOR Hamiltonian from cached Spinach operators. Syntax:
%
%      H=kehl_free_ham(parameters,paramsENDOR,spin_system,v_off_S,...
%                      spin_map,use_dipolar,term_map,euler_angles)
%
% Parameters:
%
%   parameters       - Kehl ENDOR context parameter structure.
%   paramsENDOR      - map containing ENDOR parameters.
%   spin_system      - reduced Spinach spin system.
%   v_off_S          - electron offset frequency.
%   spin_map         - ENDOR data indices represented by the reduced system.
%   use_dipolar      - flag enabling the legacy first-nucleus dipolar term.
%   term_map         - optional interaction-index override map.
%   euler_angles     - Kehl orientation angles stored by the context.
%
% Outputs:
%
%   H                - static Hilbert-space Hamiltonian.
%
% February 2024 A. Kehl (akehl@gwdg.de)
% May 2026 Spinach integration
%
% <https://spindynamics.org/wiki/index.php?title=kehl_free_ham.m>

function H=kehl_free_ham(parameters,paramsENDOR,spin_system,v_off_S,...
        spin_map,use_dipolar,term_map,euler_angles)

    % Default is to preserve the legacy dipolar correction
    if nargin<6
        use_dipolar=true;
    end

    % Default is to use the same data index for every interaction type
    if nargin<7
        term_map=struct();
    end

    % Default orientation is the input tensor frame
    if nargin<8
        euler_angles=[0 0 0];
    end

    % Check consistency
    grumble(parameters,paramsENDOR,spin_system,v_off_S,spin_map,...
        use_dipolar,term_map,euler_angles);

    % Retrieve cached Cartesian operators
    ham=kehl_ham_basis(parameters,paramsENDOR,spin_system,spin_map,...
        use_dipolar,term_map);

    % Project tensor data onto the selected Kehl orientation
    terms=kehl_orient_terms(parameters,euler_angles);

    % Get the ENDOR Larmor frequencies
    v_L=paramsENDOR('v_L');

    % Start with the electron frequency offset
    H=2*pi*v_off_S*ham.Sz;

    % Add single-nucleus terms in the legacy Kehl form
    for n=1:numel(ham.spin_map)
        zeeman_idx=ham.term_map.zeeman(n);
        hf_idx=ham.term_map.hyperfine(n);
        nqi_idx=ham.term_map.nqi(n);
        cs_idx=ham.term_map.cs(n);
        cs_larm_idx=ham.term_map.cs_larmor(n);
        H=H-2*pi*v_L(zeeman_idx)*ham.Iz{n}+...
            2*pi*v_L(cs_larm_idx)*terms.cs_zz(cs_idx)*ham.Iz{n};
        H=H+2*pi*terms.hf_zz(hf_idx)*(ham.Sz*ham.Iz{n})+...
            2*pi*terms.hf_zy(hf_idx)*(ham.Sz*ham.Iy{n})+...
            2*pi*terms.hf_zx(hf_idx)*(ham.Sz*ham.Ix{n});
        if parameters.nqi_active
            if parameters.Bterm
                H=H+nqi_full(terms.nqi(nqi_idx,:,:),ham,n);
            else
                H=H+pi*terms.nqi_zz(nqi_idx)*...
                    (3*ham.Iz{n}*ham.Iz{n}-...
                    ham.spin_q(n)*(ham.spin_q(n)+1)*ham.id);
            end
        end
    end

    % Add the legacy first-nucleus dipolar correction
    if ~isempty(ham.dip_oper)
        for n=1:numel(ham.dip_oper)
            H=H+2*pi*terms.dip_zz(n)*ham.dip_oper{n};
        end
    end

    % Return a Hermitian dense matrix for the sequence kernels
    H=full((H+H')/2);

end

% Full quadrupolar B-term contribution
function H=nqi_full(nqi_tensor,ham,spin_idx)

    % Reshape the selected tensor slice
    Q=squeeze(nqi_tensor);

    % Assemble the Cartesian quadratic form
    H=Q(1,1)*ham.Ix{spin_idx}*ham.Ix{spin_idx}+...
        Q(1,2)*ham.Ix{spin_idx}*ham.Iy{spin_idx}+...
        Q(1,3)*ham.Ix{spin_idx}*ham.Iz{spin_idx}+...
        Q(2,1)*ham.Iy{spin_idx}*ham.Ix{spin_idx}+...
        Q(2,2)*ham.Iy{spin_idx}*ham.Iy{spin_idx}+...
        Q(2,3)*ham.Iy{spin_idx}*ham.Iz{spin_idx}+...
        Q(3,1)*ham.Iz{spin_idx}*ham.Ix{spin_idx}+...
        Q(3,2)*ham.Iz{spin_idx}*ham.Iy{spin_idx}+...
        Q(3,3)*ham.Iz{spin_idx}*ham.Iz{spin_idx};

end

% Consistency enforcement
function grumble(parameters,paramsENDOR,spin_system,v_off_S,spin_map,...
        use_dipolar,term_map,euler_angles)
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
    if (~isnumeric(v_off_S))||(~isscalar(v_off_S))||(~isreal(v_off_S))
        error('v_off_S must be a real scalar.');
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
    if (~isnumeric(euler_angles))||(~isreal(euler_angles))||...
            (numel(euler_angles)~=3)
        error('euler_angles must be a three-element real vector.');
    end
end

