%KEHL_FREE_HAM Static ENDOR Hamiltonian from Spinach interactions.
%
%   H=KEHL_FREE_HAM(PARAMETERS,PARAMSENDOR,SPIN_SYSTEM,V_OFF_S,
%   SPIN_MAP,HF_ZZ,HF_ZY,HF_ZX,NQI,NQI_ZZ,CS_ZZ,D_ZZ,USE_DIPOLAR,
%   TERM_MAP) builds the orientation-selected static Kehl ENDOR
%   Hamiltonian in Hilbert space. The interaction tensors are written
%   into a temporary Spinach spin-system object and assembled with
%   hamiltonian.m and orientation.m.
%
%   Inputs:
%
%      PARAMETERS   - Kehl ENDOR parameter structure.
%
%      PARAMSENDOR  - Kehl ENDOR derived-parameter map.
%
%      SPIN_SYSTEM  - reduced Spinach spin system used by the sequence.
%
%      V_OFF_S      - electron offset frequency, Hz.
%
%      SPIN_MAP     - ENDOR data indices represented by the reduced
%                     nuclear spins in SPIN_SYSTEM.
%
%      HF_ZZ        - selected SzIz hyperfine coefficients, Hz.
%
%      HF_ZY        - selected SzIy hyperfine coefficients, Hz.
%
%      HF_ZX        - selected SzIx hyperfine coefficients, Hz.
%
%      NQI          - selected full NQI matrices, rad/s.
%
%      NQI_ZZ       - selected zz NQI coefficients.
%
%      CS_ZZ        - selected chemical-shift correction factors.
%
%      D_ZZ         - selected nuclear dipolar zz coefficients, Hz.
%
%      USE_DIPOLAR  - true to include the legacy first-nucleus dipolar
%                     correction.
%
%      TERM_MAP     - optional index overrides for legacy branches.
%
%   Output:
%
%      H            - static Hilbert-space Hamiltonian matrix.
%
%   Spinach architecture migration May 2026 Talos

function H=kehl_free_ham(parameters,paramsENDOR,spin_system,v_off_S,...
                         spin_map,HF_zz,HF_zy,HF_zx,NQI,NQI_zz,...
                         CS_zz,D_zz,use_dipolar,term_map)

% Default is to preserve the legacy dipolar correction
if nargin<13
    use_dipolar=true;
end

% Default is to use the same data index for every interaction type
if nargin<14
    term_map=struct();
end

% Check consistency
grumble(parameters,paramsENDOR,spin_system,v_off_S,spin_map,HF_zz,...
        HF_zy,HF_zx,NQI,NQI_zz,CS_zz,D_zz,use_dipolar,term_map);

% Get the ENDOR Larmor frequencies
v_L=paramsENDOR('v_L');

% Normalise the spin-map shape
spin_map=spin_map(:).';
n_nuc=numel(spin_map);
n_spins=spin_system.comp.nspins;

% Expand legacy index overrides
term_map=complete_map(term_map,spin_map);

% Reset Zeeman tensors in the temporary spin system
spin_system.inter.zeeman.matrix=cell(1,n_spins);
for n=1:n_spins
    spin_system.inter.zeeman.matrix{n}=zeros(3,3);
end

% Reset coupling tensors in the temporary spin system
spin_system.inter.coupling.matrix=cell(n_spins,n_spins);
for n=1:n_spins
    for k=1:n_spins
        spin_system.inter.coupling.matrix{n,k}=[];
    end
end

% Install the electron offset term
spin_system.inter.zeeman.matrix{1}=eye(3,3)*(2*pi*v_off_S);

% Install selected nuclear terms
for n=1:n_nuc
    spin_idx=n+1;
    zeeman_idx=term_map.zeeman(n);
    hf_idx=term_map.hyperfine(n);
    nqi_idx=term_map.nqi(n);
    cs_idx=term_map.cs(n);
    cs_larmor_idx=term_map.cs_larmor(n);

    % Install the nuclear Zeeman and chemical-shift contribution
    spin_system.inter.zeeman.matrix{spin_idx}=...
        eye(3,3)*(-2*pi*v_L(zeeman_idx)+...
                  2*pi*v_L(cs_larmor_idx)*CS_zz(cs_idx));

    % Install the left-secular hyperfine tensor
    spin_system.inter.coupling.matrix{1,spin_idx}=...
        2*pi*[0 0 HF_zx(hf_idx);...
              0 0 HF_zy(hf_idx);...
              0 0 HF_zz(hf_idx)];

    % Install either full or secularised NQI terms
    if parameters.Bterm==true
        spin_system.inter.coupling.matrix{spin_idx,spin_idx}=...
            squeeze(NQI(nqi_idx,:,:));
    else
        spin_system.inter.coupling.matrix{spin_idx,spin_idx}=...
            diag([-pi*NQI_zz(nqi_idx) -pi*NQI_zz(nqi_idx) ...
                   2*pi*NQI_zz(nqi_idx)]);
    end
end

% Install the legacy first-nucleus dipolar correction
if use_dipolar&&parameters.dipolar_active&&(n_nuc>1)
    for n=2:min(n_nuc,numel(D_zz)+1)
        spin_idx=n+1;
        dip_coeff=2*pi*D_zz(n-1);
        dip_matrix=-dip_coeff*ones(3,3)/2;
        dip_matrix(3,3)=dip_matrix(3,3)+3*dip_coeff/2;
        spin_system.inter.coupling.matrix{2,spin_idx}=dip_matrix;
    end
end

% Set laboratory-frame assumptions for the synthetic tensors
spin_system=assume(spin_system,'labframe');

% Keep only electron-z terms from the hyperfine tensors
for n=1:n_nuc
    spin_system=dictum(spin_system,[1 n+1],'z*');
end

% Assemble the full static Hamiltonian through Spinach
[H_iso,H_aniso]=hamiltonian(spin_system);
H=full(H_iso+orientation(H_aniso,[0 0 0]));

end

% Consistency enforcement
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
function grumble(parameters,paramsENDOR,spin_system,v_off_S,spin_map,...
                 HF_zz,HF_zy,HF_zx,NQI,NQI_zz,CS_zz,D_zz,...
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
if (~isnumeric(HF_zz))||(~isnumeric(HF_zy))||(~isnumeric(HF_zx))
    error('hyperfine coefficient arrays must be numeric.');
end
if (~isnumeric(NQI))||(~isnumeric(NQI_zz))||(~isnumeric(CS_zz))
    error('NQI and chemical-shift coefficient arrays must be numeric.');
end
if ~isnumeric(D_zz)
    error('D_zz must be numeric.');
end
if (~islogical(use_dipolar))&&(~isnumeric(use_dipolar))
    error('use_dipolar must be logical or numeric.');
end
if ~isstruct(term_map)
    error('term_map must be a structure.');
end
end
