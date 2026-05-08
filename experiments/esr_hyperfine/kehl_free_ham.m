%KEHL_FREE_HAM Static ENDOR Hamiltonian from cached Spinach basis.
%
%   H=KEHL_FREE_HAM(PARAMETERS,PARAMSENDOR,SPIN_SYSTEM,V_OFF_S,
%   SPIN_MAP,USE_DIPOLAR,TERM_MAP,EULER_ANGLES) assembles the
%   orientation-specific static Kehl ENDOR Hamiltonian in Hilbert space.
%   The expensive isotropic Hamiltonian and rotational basis are cached by
%   kehl_ham_basis.m; this function only combines them with orientation.m
%   and adds the electron offset term.
%
%   Inputs:
%
%      PARAMETERS    - Kehl ENDOR parameter structure.
%
%      PARAMSENDOR   - Kehl ENDOR derived-parameter map.
%
%      SPIN_SYSTEM   - reduced Spinach spin system used by the sequence.
%
%      V_OFF_S       - electron offset frequency, Hz.
%
%      SPIN_MAP      - ENDOR data indices represented by the reduced
%                      nuclear spins in SPIN_SYSTEM.
%
%      USE_DIPOLAR   - true to include the legacy first-nucleus dipolar
%                      correction.
%
%      TERM_MAP      - optional index overrides for legacy branches.
%
%      EULER_ANGLES  - Spinach ZYZ orientation angles, radians.
%
%   Output:
%
%      H             - static Hilbert-space Hamiltonian matrix.
%
%   Spinach architecture migration May 2026 Talos

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

% Retrieve the cached Spinach Hamiltonian basis
ham=kehl_ham_basis(parameters,paramsENDOR,spin_system,spin_map,...
                   use_dipolar,term_map);

% Assemble the orientation-specific Hamiltonian
H=ham.I+orientation(ham.Q,euler_angles)+2*pi*v_off_S*ham.Sz;

% Add legacy dipolar corrections by cheap per-orientation trigonometry
if ~isempty(ham.dip_oper)
    R=euler2dcm(euler_angles);
    for n=1:numel(ham.dip_oper)
        D=R*ham.dip_matrix{n}*R';
        H=H+2*pi*D(3,3)*ham.dip_oper{n};
    end
end

% Return a Hermitian dense matrix for the sequence kernels
H=full((H+H')/2);

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
