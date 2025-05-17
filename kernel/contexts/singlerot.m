% Single angle spinning context. In Liouville space, this wrapper builds
% the Fokker-Planck evolution generator that includes the spin Hamiltoni-
% an commutation superoperator, applicable dissipators (relaxation, kine-
% tics), and the rotor turning generator. In Hilbert space, this wrapper
% builds the stack of spin Hamiltonians, one for each rotor phase. Those
% are handed over to the pulse sequence, which the user must supply as a
% function handle. Syntax:
%
%       [answer,sph_grid]=singlerot(spin_system,pulse_sequence,...
%                                   parameters,assumptions)
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
%                         3-element vector; this is the direction
%                         around which the rotor is turning
%
%   parameters.spins    - a cell array giving the spins that 
%                         the pulse sequence involves, e.g. 
%                         {'1H','13C'}
%
%   parameters.offset   - a cell array giving transmitter off-
%                         sets in Hz on each of the spins listed
%                         in parameters.spins array
%
%   parameters.max_rank - maximum rotor harmonic rank to retain
%                         in the solution (increase till conver-
%                         gence is achieved, a good guess value
%                         is the number of spinning sidebands
%                         expected in the spectrum)
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
%   parameters.grid     - spherical grid file name; see grids
%                         directory in the kernel. Two-angle
%                         grids should be used in Liouville
%                         space and three-angle grids in Hil-
%                         bert space.
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
% Outputs:
%
%     answer - the poweder average or a cell array ofwhatever it is 
%              that the pulse sequence returns
%
%     sph_grid - spherical grid used ithe calculation
%
% Note: arbitrary order rotating frame transformation is supported, inc-
%       luding infinite order. See the header of rotframe.m for further
%       information.
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=singlerot.m>

function [answer,sph_grid]=singlerot(spin_system,pulse_sequence,...
                                     parameters,assumptions)

% Show the banner
banner(spin_system,'sequence_banner'); 

% Set common defaults
parameters=defaults(spin_system,parameters);

% Check consistency
grumble(spin_system,pulse_sequence,parameters,assumptions);

% Load the spherical integration grid
sph_grid=load([spin_system.sys.root_dir '/kernel/grids/' ...
               parameters.grid],'alphas','betas','gammas','weights');
alphas=sph_grid.alphas; betas=sph_grid.betas; 
gammas=sph_grid.gammas; weights=sph_grid.weights;

% Get and report rotor axis orientation angles
[rotor_phi,rotor_theta,~]=cart2sph(parameters.axis(1),...
                                   parameters.axis(2),...
                                   parameters.axis(3)); 
rotor_theta=pi/2-rotor_theta;
report(spin_system,['lab frame rotor direction, theta: ' num2str(180*rotor_theta/pi) ' degrees.']);
report(spin_system,['lab frame rotor direction, phi:   ' num2str(180*rotor_phi/pi) ' degrees.']);

% Report problem dimensions to the user
spn_dim=size(spin_system.bas.basis,1);
spc_dim=2*parameters.max_rank+1; n_orients=numel(weights);
report(spin_system,['spin space problem dimension:     ' num2str(spn_dim)]);
report(spin_system,['number of rotor grid points:      ' num2str(spc_dim)]);
report(spin_system,['number of powder grid points:     ' num2str(n_orients)]);

% Forward dimension information to the pulse sequence
parameters.spc_dim=spc_dim; parameters.spn_dim=spn_dim;

% Set the interaction assumptions
spin_system=assume(spin_system,assumptions);

% Get isotropic Hamiltonian and rotations basis
report(spin_system,'building the Hamiltonian...');
[I,Q]=hamiltonian(assume(spin_system,assumptions));

% Apply channel frequency offsets
report(spin_system,'applying frequency offsets...');
I=frqoffset(spin_system,I,parameters);

% Compute isotropic thermal equilibrium
if ismember('iso_eq',parameters.needs)
    report(spin_system,'WARNING - equilibrium state uses the isotropic Hamitonian.');
    I_labframe=hamiltonian(assume(spin_system,'labframe'),'left');
    parameters.rho0=equilibrium(spin_system,I_labframe);
end

% Get carrier operators for numerical
% rotating frame transformations
C=cell(size(parameters.rframes));
for n=1:numel(parameters.rframes)
    C{n}=carrier(spin_system,parameters.rframes{n}{1});
end

% Get relaxation and kinetics generators
R=relaxation(spin_system); K=kinetics(spin_system);

% Formalism-dependent stage
switch spin_system.bas.formalism

    % Liouville space
    case {'sphten-liouv','zeeman-liouv'}

        % Report the composite dimension of the Fokker-Planck problem
        report(spin_system,['Fokker-Planck problem dimension:  ' num2str(spc_dim*spn_dim)]);

        % Make the rotor turning generator
        [rotor_phases,d_dphi]=fourdif(spc_dim,1);
        M=2*pi*parameters.rate*kron(d_dphi,speye([spn_dim spn_dim]));

        % Project relaxation and kinetics superoperators into the FP space
        R=kron(speye([spc_dim spc_dim]),R); K=kron(speye([spc_dim spc_dim]),K);

        % Project the initial state into the FP space
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

        % Project the coil state
        if isfield(parameters,'coil')
            space_part=ones(parameters.spc_dim,1);
            parameters.coil=kron(space_part,parameters.coil);
        end

    % Hilbert space
    case {'zeeman-hilb','zeeman-wavef'}

        % Get rotor phases and avoid parfor bug
        rotor_phases=fourdif(spc_dim,1); M=[];

    otherwise

        % Complain and bomb out
        error('unknown formalism specification.');

end

% Preallocate answer array
ans_array=cell(numel(weights),1);

% Decide the parallelisation strategy
if isfield(parameters,'serial')&&parameters.serial
       
    % Serial execution at this level with a warning to the console
    nworkers=0; report(spin_system,'WARNING: powder grid parallelisation is turned off.');
    
else
    
    % Parallel execution
    nworkers=min([poolsize n_orients]);

end

% Inform the user and silence the output
prev_setting=spin_system.sys.output; ss_prev.sys.output=spin_system.sys.output;
if ~isfield(parameters,'verbose')||(parameters.verbose==0)
    report(spin_system,'pulse sequence silenced to avoid excessive output.')
    spin_system.sys.output='hush';
end

% Parfor rigging
if ~isworkernode
    DQ=parallel.pool.DataQueue;
    afterEach(DQ,@(~)parfor_progr);
    orients_done=0; last_toc=0;
    tic; ticBytes(gcp); do_diag=true;
else
    do_diag=false; DQ=[];
end

% Parfor progress updater
function parfor_progr()
    orients_done=orients_done+1; last_message=toc-last_toc;
    if (last_message>5)||(orients_done==n_orients)
        report(ss_prev,[num2str(orients_done) '/' num2str(n_orients) ' orientations done, ' ...
                        num2str(orients_done/toc) ' orientations per second.']); 
        last_toc=toc;
    end
end

% Parallel powder averaging loop
parfor (q=1:n_orients,nworkers) %#ok<*PFBNS>

    % Preallocate Hamiltonian blocks
    H=cell(2*parameters.max_rank+1,1); H(:)={I};

    % Build Hamiltonian rotor stack
    for n=1:(2*parameters.max_rank+1)
        
        % Loop over spherical ranks
        for r=1:numel(Q)
            
            % Compute crystallite orientation
            D_mol2rot=wigner(r,alphas(q),betas(q),gammas(q));
            
            % Compute rotor axis tilt
            D_lab2rot=wigner(r,rotor_phi,rotor_theta,0);
        
            % Compute rotor rotation
            D_rotor=wigner(r,0,0,rotor_phases(n));
            
            % Compose rotations
            D_comp=D_lab2rot*D_rotor*D_mol2rot;
        
            % Build the block
            for k=1:(2*r+1)
                for m=1:(2*r+1)
                    H{n}=H{n}+D_comp(k,m)*Q{r}{k,m};
                end
            end
            
        end
        
        % Apply rotating frames
        for k=1:numel(parameters.rframes)

            % Arithmetic clean-up
            H{n}=(H{n}+H{n}')/2;
            
            % Rotating frame transformation
            H{n}=rotframe(spin_system,C{k},H{n},...
                          parameters.rframes{k}{1},...
                          parameters.rframes{k}{2});

        end

        % H must be sparse
        H{n}=sparse(H{n});
        
    end
    
    % Formalism-dependent stage
    switch spin_system.bas.formalism

        % Liouville space
        case {'sphten-liouv','zeeman-liouv'}

            % Assemble the Fokker-Planck evolution generator
            G=clean_up(spin_system,blkdiag(H{:})+1i*M,spin_system.tols.liouv_zero);
    
            % Run the pulse sequence
            ans_array{q}=pulse_sequence(spin_system,parameters,G,R,K);

        % Hilbert space
        case {'zeeman-hilb','zeeman-wavef'}

            % Run the pulse sequence with a Hamiltonian stack
            ans_array{q}=pulse_sequence(spin_system,parameters,H,R,K);

        otherwise

            % Complain and bomb out
            error('unknown formalism specification.');

    end

    % Report parfor progress
    if do_diag, send(DQ,n); end
    
end

% Unsilence the output
spin_system.sys.output=prev_setting;

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
    
    % Return components
    answer=ans_array;
    
    % Inform the user
    report(spin_system,'returning pulse sequence outputs at each orientation...');
    
end

% Parfor performance
if ~isworkernode
    nbytes=mean(tocBytes(gcp),1)/2^20; walltime=toc();
    report(spin_system,['average worker process received ' num2str(nbytes(1)) ...
                        ' MB and sent back ' num2str(nbytes(2)) ' MB']);
    report(spin_system,['powder average run time: ' num2str(walltime) ' seconds']);
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

% Option combination restrictions
if isfield(parameters,'rho0')&&isfield(parameters,'needs')&&...
   ismember('iso_eq',parameters.needs)
    error('thermal equilibrium request conflicts with user-specified initial condition.');
end
if strcmp(spin_system.bas.formalism,'zeeman-wavef')&&isfield(parameters,'needs')&&...
   ismember('iso_eq',parameters.needs)
    error('thermal equilibrium state cannot be represented by a wavefunction.');
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

