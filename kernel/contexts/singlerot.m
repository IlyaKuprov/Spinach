% Single angle spinning context. Generates a Liouvillian superoperator 
% and passes it on to the pulse sequence function, which should be sup-
% plied as a handle. Syntax:
%
%  answer=singlerot(spin_system,pulse_sequence,parameters,assumptions)
%
% Parameters:
%
%   pulse_sequence      - pulse sequence function handle. See the
%                         experiments directory for the list of
%                         pulse sequences that ship with Spinach.
%
%   parameters.rate     - spinning rate in Hz. Positive numbers
%                         for JEOL, negative for Varian and Bruker
%                         due to different rotation directions.
%
%   parameters.axis     - spinning axis, given as a normalized
%                         3-element vector
%
%   parameters.spins    - a cell array giving the spins that 
%                         the pulse sequence involves, e.g. 
%                         {'1H','13C'}
%
%   parameters.offset   - a cell array giving transmitter off-
%                         sets in Hz on each of the spins listed
%                         in parameters.spins array
%
%   parameters.max_rank - maximum harmonic rank to retain in
%                         the solution (increase till conver-
%                         gence is achieved, approximately
%                         equal to the number of spinning si-
%                         debands in the spectrum)
%
%   parameters.rframes  - rotating frame specification, e.g.
%                         {{'13C',2},{'14N',3}} requests second
%                         order rotating frame transformation
%                         with respect to carbon-13 and third
%                         order rotating frame transformation
%                         with respect to nitrogen-14. When
%                         this option is used, the assumptions
%                         on the respective spins should be
%                         laboratory frame.
%
%   parameters.grid     - spherical grid file name. See grids
%                         directory in the kernel.
%
%   parameters.needs    - a cell array of character strings spe-
%                         cifying additional requirements that
%                         the sequence has:
%
%                          'iso_eq' - thermal equilibrium state
%                          of the isotropic Hamiltonian will be
%                          placed into parameters.rho0
%
%   parameters.sum_up   - when set to 1 (default), returns the
%                         powder average. When set to 0, returns
%                         individual answers for each point in 
%                         the powder as a cell array.
%
%   parameters.*        - additional subfields may be required 
%                         by the pulse sequence - check its do-
%                         cumentation page 
%
%   assumptions         - context-specific assumptions ('nmr', 'epr',
%                         'labframe', etc.) - see the pulse sequence
%                         header for information on this setting.
%
% The parameters structure is passed to the pulse sequence with the follo-
% wing additional parameters set:
%
%   parameters.spc_dim  - matrix dimension for the spatial
%                         dynamics subspace
%
%   parameters.spn_dim  - matrix dimension for the spin 
%                         dynamics subspace
%
% This function returns the powder average of whatever it is that the pulse
% sequence returns.
%
% Note: arbitrary order rotating frame transformation is supported, inc-
%       luding infinite order. See the header of rotframe.m for further
%       information.
%
% Note: the function supports parallel processing via Matlab's Distri-
%       buted Computing Toolbox - different system orientations are eva-
%       luated on different Matlab workers.
%
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=singlerot.m>

function answer=singlerot(spin_system,pulse_sequence,parameters,assumptions)

% Show the banner
banner(spin_system,'sequence_banner'); 

% Set common defaults
parameters=defaults(spin_system,parameters);

% Check consistency
grumble(spin_system,pulse_sequence,parameters,assumptions);

% Report to the user
report(spin_system,'building the Liouvillian...');

% Set the assumptions
spin_system=assume(spin_system,assumptions);

% Get the Hamiltonian
[H,Q]=hamiltonian(spin_system);

% Apply offsets
H=frqoffset(spin_system,H,parameters);

% Compute the thermal equilibrium
if ismember('iso_eq',parameters.needs)
    report(spin_system,'WARNING - thermal equilibrium uses the isotropic Hamitonian.');
    if isfield(parameters,'rho0')
        report(spin_system,'WARNING - user-specified initial condition has been ignored.');
    end
    parameters.rho0=equilibrium(spin_system,hamiltonian(assume(spin_system,'labframe'),'left'));
end

% Get problem dimensions
spc_dim=2*parameters.max_rank+1; spn_dim=size(H,1);
report(spin_system,['lab space problem dimension     ' num2str(spc_dim)]);
report(spin_system,['spin space problem dimension    ' num2str(spn_dim)]);
report(spin_system,['Fokker-Planck problem dimension ' num2str(spc_dim*spn_dim)]);
parameters.spc_dim=spc_dim; parameters.spn_dim=spn_dim;

% Compute rotor angles and the derivative operator
[rotor_angles,d_dphi]=fourdif(spc_dim,1);

% Compute spinning operator
M=2*pi*parameters.rate*kron(d_dphi,speye(size(H)));

% Get rotor axis orientation
[rotor_phi,rotor_theta,~]=cart2sph(parameters.axis(1),...
                                   parameters.axis(2),...
                                   parameters.axis(3)); 
rotor_theta=pi/2-rotor_theta; 

% Get carrier operators
C=cell(size(parameters.rframes));
for n=1:numel(parameters.rframes)
    C{n}=carrier(spin_system,parameters.rframes{n}{1});
end

% Get relaxation and kinetics 
R=relaxation(spin_system); K=kinetics(spin_system);

% Get the averaging grid as a structure
sph_grid=load([spin_system.sys.root_dir '/kernel/grids/' ...
               parameters.grid],'alphas','betas','gammas','weights');
           
% Assign local variables
alphas=sph_grid.alphas; betas=sph_grid.betas; 
gammas=sph_grid.gammas; weights=sph_grid.weights;

% Project relaxation and kinetics
R=kron(speye(spc_dim),R); K=kron(speye(spc_dim),K);

% Project the initial state
if isfield(parameters,'rho0')&&strcmp(parameters.grid,'single_crystal')
    
    % Single crystal simulations start at 12 o'clock
    space_part=zeros(parameters.spc_dim,1); space_part(1)=1;
    parameters.rho0=kron(space_part,parameters.rho0);
    report(spin_system,'single crystal simulation, rotor phase averaging switched off.');
    
elseif isfield(parameters,'rho0')
    
    % Powder simulations start equally distributed
    space_part=ones(parameters.spc_dim,1)/parameters.spc_dim;
    parameters.rho0=kron(space_part,parameters.rho0);
    report(spin_system,'powder simulation, rotor phase averaging switched on.');
    
end

% Project the coil state: same coil everywhere
if isfield(parameters,'coil')
    space_part=ones(parameters.spc_dim,1);
    parameters.coil=kron(space_part,parameters.coil);
end

% Preallocate answer array
ans_array=cell(numel(weights),1);

% Run serially if needed
if isfield(parameters,'serial')&&...
           parameters.serial
       
    % Serial execution
    nworkers=0;
    
    % Inform the user
    report(spin_system,'WARNING: parallelisation turned off by the user.');
    
else
    
    % Parallel execution
    nworkers=min([poolsize numel(weights)]);

end

% Inform the user and silence the output
prev_setting=spin_system.sys.output;
report(spin_system,['powder simulation with ' num2str(numel(weights)) ' orientations.']);
if ~isfield(parameters,'verbose')||(parameters.verbose==0)
    report(spin_system,'pulse sequence silenced to avoid excessive output.')
    spin_system.sys.output='hush';
end

% MDCS diagnostics
parallel_profiler_start;

% Powder averaging loop
parfor (q=1:numel(weights),nworkers) %#ok<*PFBNS>

    % Preallocate Liouvillian blocks
    L=cell(2*parameters.max_rank+1,...
           2*parameters.max_rank+1); 
    for n=1:(2*parameters.max_rank+1)
        for k=1:(2*parameters.max_rank+1)
            L{n,k}=krondelta(n,k)*H;
        end
    end

    % Build Liouvillian blocks
    for n=1:(2*parameters.max_rank+1)
        
        % Loop over spherical ranks
        for r=1:numel(Q)
            
            % Compute crystallite orientation
            D_mol2rot=wigner(r,alphas(q),betas(q),gammas(q));
            
            % Compute rotor axis tilt
            D_lab2rot=wigner(r,rotor_phi,rotor_theta,0);
        
            % Compute rotor rotation
            D_rotor=wigner(r,0,0,rotor_angles(n));
            
            % Compose rotations
            D=D_lab2rot*D_rotor*D_mol2rot;
        
            % Build the block
            for k=1:(2*r+1)
                for m=1:(2*r+1)
                    L{n,n}=L{n,n}+D(k,m)*Q{r}{k,m};
                end
            end
            
        end
        
        % Apply rotating frames
        for k=1:numel(parameters.rframes)
            L{n,n}=rotframe(spin_system,C{k},(L{n,n}+L{n,n}')/2,...
                            parameters.rframes{k}{1},parameters.rframes{k}{2});
        end
        
    end

    % Assemble the Liouvillian
    L=clean_up(spin_system,cell2mat(L)+1i*M,spin_system.tols.liouv_zero);
    
    % Run the pulse sequence
    ans_array{q}=pulse_sequence(spin_system,parameters,L,R,K);
    
end

% Unsilence the output
spin_system.sys.output=prev_setting;

% Get MDCS diagnostics
parallel_profiler_report;

% Decide the return array
if parameters.sum_up
    
    % Return weighted sum
    answer=weights(1)*ans_array{1};
    for n=2:numel(ans_array)
        answer=answer+weights(n)*ans_array{n};
    end
    
    % Inform the user
    report(spin_system,'returning powder averaged pulse sequence output...');
    
else
    
    % Return components and weights
    answer.components=ans_array;
    answer.weights=weights;
    
    % Inform the user
    report(spin_system,'returning pulse sequence outputs at each orientation...');
    
end

end

% Default parameters
function parameters=defaults(spin_system,parameters)
if ~isfield(parameters,'decouple')
    report(spin_system,'parameters.decouple field not set, assuming no decoupling.');
    parameters.decouple={};
end
if ~isfield(parameters,'rframes')
    report(spin_system,'parameters.rframes field not set, assuming no additional rotating frame transformations.');
    parameters.rframes={};
end
if ~isfield(parameters,'offset')
    report(spin_system,'parameters.offset field not set, assuming zero offsets.');
    parameters.offset=zeros(size(parameters.spins));
end
if ~isfield(parameters,'verbose')
    report(spin_system,'parameters.verbose field not set, silencing array operations.');
    parameters.verbose=0;
end
if ~isfield(parameters,'sum_up')
    parameters.sum_up=1;
end
if ~isfield(parameters,'needs')
    parameters.needs={};
end
end

% Consistency enforcement
function grumble(spin_system,pulse_sequence,parameters,assumptions)

% Formalism 
if ~ismember(spin_system.bas.formalism,{'zeeman-liouv','sphten-liouv'})
    error('this function is only available in Liouville space.');
end

% Rotor rank
if ~isfield(parameters,'max_rank')
    error('parameters.max_rank subfield must be present.');
elseif (~isnumeric(parameters.max_rank))||(~isreal(parameters.max_rank))||...
       (mod(parameters.max_rank,1)~=0)||(parameters.max_rank<0)
    error('parameters.max_rank must be a positive real integer.');
end

% Spinning rate
if ~isfield(parameters,'rate')
    error('parameters.rate subfield must be present.');
elseif (~isnumeric(parameters.rate))||(~isreal(parameters.rate))
    error('parameters.rate must be a real number.');
end

% Spinning axis
if ~isfield(parameters,'axis')
    error('parameters.axis subfield must be present.');
elseif (~isnumeric(parameters.axis))||(~isreal(parameters.axis))||...
       (~isrow(parameters.axis))||(numel(parameters.axis)~=3)
    error('parameters.axis must be a row vector of three real numbers.');
end

% Spherical grid
if ~isfield(parameters,'grid')
    error('spherical averaging grid must be specified in parameters.grid variable.');
elseif isempty(parameters.grid)
    error('parameters.grid variable cannot be empty.');
elseif ~ischar(parameters.grid)
    error('parameters.grid variable must be a character string.');
end

% Pulse sequence
if ~isa(pulse_sequence,'function_handle')
    error('pulse_sequence argument must be a function handle.');
end

% Assumptions
if ~ischar(assumptions)
    error('assumptions argument must be a character string.');
end

% Active spins
if isempty(parameters.spins)
    error('parameters.spins variable cannot be empty.');
elseif ~iscell(parameters.spins)
    error('parameters.spins variable must be a cell array.');
elseif ~all(cellfun(@ischar,parameters.spins))
    error('all elements of parameters.spins cell array must be strings.');
elseif any(~ismember(parameters.spins,spin_system.comp.isotopes))
    error('parameters.spins refers to a spin that is not present in the system.');
end

% Offsets
if isempty(parameters.offset)
    error('parameters.offset variable cannot be empty.');
elseif ~isnumeric(parameters.offset)
    error('parameters.offset variable must be an array of real numbers.');
elseif ~isfield(parameters,'spins')
    error('parameters.spins variable must be specified alongside parameters.offset variable.');
elseif numel(parameters.offset)~=numel(parameters.spins)
    error('parameters.offset variable must have the same number of elements as parameters.spins.');
end

% Rotating frame transformations
if ~isfield(parameters,'rframes')
    error('parameters.rframes variable must be specified.');
elseif ~iscell(parameters.rframes)
    error('parameters.rframes must be a cell array.');
end
for n=1:numel(parameters.rframes)
    if ~iscell(parameters.rframes{n})
        error('elements of parameters.rframes must be cell arrays.');
    end
    if numel(parameters.rframes{n})~=2
        error('elements of parameters.rframes must have exactly two sub-elements each.');
    end
    if ~ischar(parameters.rframes{n}{1})
        error('the first part of each element of parameters.rframes must be a character string.');
    end
    if ~ismember(parameters.rframes{n}{1},spin_system.comp.isotopes)
        error('parameters.rframes refers to a spin that is not present in the system.');
    end
    if ~isnumeric(parameters.rframes{n}{2})
        error('the second part of each element of parameters.rframes must be a number.');
    end
end
        
end

% There are no set working hours, your contract requires you to work such
% hours as are reasonably necessary to perform your duties.
%
% [...]
%
% Associate Professors are not permitted to use the designation Professor.
% Any Associate Professor who mistakenly uses the title Professor will be
% told that this is not correct and asked to amend it.
% 
% IK's contract at Southampton University, 2014

