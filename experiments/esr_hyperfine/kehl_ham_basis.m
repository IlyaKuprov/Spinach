%KEHL_HAM_BASIS Cached Spinach Hamiltonian basis for Kehl kernels.
%
%   HAM=KEHL_HAM_BASIS(PARAMETERS,PARAMSENDOR,SPIN_SYSTEM,SPIN_MAP,
%   USE_DIPOLAR,TERM_MAP) returns a cached Hamiltonian basis for the
%   reduced Kehl ENDOR spin system. The isotropic part HAM.I and the
%   rotational basis HAM.Q are built once with hamiltonian.m; orientation-
%   specific Hamiltonians are then assembled with orientation.m.
%
%   Inputs:
%
%      PARAMETERS  - Kehl ENDOR context parameter structure.
%
%      PARAMSENDOR - Kehl ENDOR parameter map.
%
%      SPIN_SYSTEM - reduced Spinach spin system used by the sequence.
%
%      SPIN_MAP    - ENDOR data indices represented in SPIN_SYSTEM.
%
%      USE_DIPOLAR - true to include the legacy first-nucleus dipolar
%                    correction basis.
%
%      TERM_MAP    - optional index overrides for legacy branches.
%
%   Output:
%
%      HAM         - structure containing I, Q, Sz, and metadata.
%
%   Spinach architecture migration May 2026 Talos

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
cache_key=ham_key(parameters,paramsENDOR,spin_system,spin_map,...
                  use_dipolar,term_map);

% Return the cached basis if it already exists
if isKey(basis_cache,cache_key)
    ham=basis_cache(cache_key);
    return
end

% Get the ENDOR Larmor frequencies
v_L=paramsENDOR('v_L');
n_nuc=numel(spin_map);
n_spins=spin_system.comp.nspins;

% Remove constructor carrier frequencies from the synthetic Hamiltonian
spin_system.inter.basefrqs=zeros(size(spin_system.inter.basefrqs));

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

% Install nuclear Zeeman and chemical-shift tensors
for n=1:n_nuc
    spin_idx=n+1;
    zeeman_idx=term_map.zeeman(n);
    cs_idx=term_map.cs(n);
    cs_larm_idx=term_map.cs_larmor(n);
    cs_matrix=zeros(3,3);
    if isfield(parameters,'cs_matrix')&&parameters.cs_active
        cs_matrix=parameters.cs_matrix(3*cs_idx-2:3*cs_idx,:)*1e-12;
    end
    spin_system.inter.zeeman.matrix{spin_idx}=...
        -2*pi*v_L(zeeman_idx)*eye(3,3)+...
         2*pi*v_L(cs_larm_idx)*cs_matrix;
end

% Install left-secular hyperfine tensors
for n=1:n_nuc
    spin_idx=n+1;
    hf_idx=term_map.hyperfine(n);
    A=parameters.hfc_matrix(3*hf_idx-2:3*hf_idx,:);
    spin_system.inter.coupling.matrix{1,spin_idx}=2*pi*A;
end

% Install nuclear quadrupolar tensors
if parameters.nqi_active
    for n=1:n_nuc
        spin_idx=n+1;
        nqi_idx=term_map.nqi(n);
        Q=parameters.nqi_matrix(3*nqi_idx-2:3*nqi_idx,:);
        spin_system.inter.coupling.matrix{spin_idx,spin_idx}=2*pi*Q;
    end
end

% Set laboratory-frame assumptions for the synthetic tensors
spin_system=assume(spin_system,'labframe');

% Keep only nuclear-z terms from nuclear Zeeman tensors
for n=1:n_nuc
    spin_system=dictum(spin_system,n+1,'secular');
end

% Keep only electron-z terms from the hyperfine tensors
for n=1:n_nuc
    spin_system=dictum(spin_system,[1 n+1],'z*');
end

% Use the secular NQI projection unless explicit B terms are requested
if parameters.nqi_active&&(~parameters.Bterm)
    for n=1:n_nuc
        spin_system=dictum(spin_system,[n+1 n+1],'secular');
    end
end

% Build and store the reusable Hamiltonian components
[ham.I,ham.Q]=hamiltonian(spin_system);
ham.Sz=full(operator(spin_system,'Lz',1));
ham.dip_matrix={};
ham.dip_oper={};
if use_dipolar&&parameters.dipolar_active&&(n_nuc>1)
    for n=2:min(n_nuc,size(parameters.dipolar_matrix,1)/3+1)
        spin_idx=n+1;
        ham.dip_matrix{n-1}=parameters.dipolar_matrix(3*n-5:3*n-3,:);
        ham.dip_oper{n-1}=legacy_dip_oper(spin_system,2,spin_idx);
    end
end
ham.spin_map=spin_map;
ham.use_dipolar=logical(use_dipolar);
ham.term_map=term_map;

% Store the basis in the process-local cache
basis_cache(cache_key)=ham;

end

% Legacy secular dipolar operator shape
function D_oper=legacy_dip_oper(spin_system,spin_a,spin_b)

% Set Cartesian operator labels
ops={'Lx','Ly','Lz'};

% Set the legacy axial matrix shape
dip_shape=-ones(3,3)/2;
dip_shape(3,3)=dip_shape(3,3)+3/2;

% Assemble the reusable operator combination
D_oper=sparse(size(mprealloc(spin_system,0),1),...
              size(mprealloc(spin_system,0),2));
for n=1:3
    for k=1:3
        D_oper=D_oper+dip_shape(n,k)*operator(spin_system,...
                           {ops{n},ops{k}},{spin_a,spin_b});
    end
end

% Return a dense operator for the small Kehl Hilbert spaces
D_oper=full(D_oper);

end

% Build a conservative cache key
function cache_key=ham_key(parameters,paramsENDOR,spin_system,spin_map,...
                           use_dipolar,term_map)

% Get optional interaction blocks
if isfield(parameters,'cs_matrix')
    cs_matrix=parameters.cs_matrix;
else
    cs_matrix=zeros(0,3);
end
if isfield(parameters,'dipolar_matrix')
    dip_matrix=parameters.dipolar_matrix;
else
    dip_matrix=zeros(0,3);
end

% Identify the reduced Spinach basis
if isfield(spin_system.bas,'basis_hash')
    basis_id=spin_system.bas.basis_hash;
else
    basis_id=spin_system.bas.basis;
end

% Hash all data that can affect the Hamiltonian basis
cache_key=md5_hash({spin_system.comp.isotopes,basis_id,...
                    paramsENDOR('v_L'),parameters.hfc_matrix,...
                    parameters.nqi_matrix,cs_matrix,dip_matrix,...
                    parameters.nqi_active,parameters.cs_active,...
                    parameters.dipolar_active,parameters.Bterm,...
                    spin_map,logical(use_dipolar),term_map});

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
