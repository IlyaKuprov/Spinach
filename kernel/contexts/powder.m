% Static powder interface to pulse sequences. Generates a Liouvillian 
% superoperator, the initial state, the coil state, then passes them 
% on to the pulse sequence function. Syntax:
%
%      [answer,sph_grid]=powder(spin_system,pulse_sequence,...
%                               parameters,assumptions)
%
% Parameters:
%
%   pulse_sequence     - pulse sequence function handle. See the
%                        experiments directory for the list of
%                        pulse sequences that ship with Spinach.
%  
%   parameters.spins   - a cell array giving the spins that the
%                        pulse sequence works on, in the order
%                        of channels, e.g. {'1H','13C'}
%
%   parameters.offset  - a cell array giving transmitter offsets
%                        in Hz on each of the spins listed in
%                        parameters.spins
% 
%   parameters.grid    - name of the spherical averaging grid 
%                        file (see the grids directory in the
%                        kernel).
% 
%   parameters.rframes - rotating frame specification, e.g.
%                        {{'13C',2},{'14N,3}} requests second
%                        order rotating frame transformation
%                        with respect to carbon-13 and third
%                        order rotating frame transformation
%                        with respect to nitrogen-14. When
%                        this option is used, the assumptions
%                        on the respective spins should be
%                        laboratory frame.
%
%   parameters.needs   - a cell array of strings specifying ad-
%                        ditional information required by the
%                        sequence:
% 
%                        'zeeman_op' - Zeeman part of the Hami-
%                        ltonian in the laboratory frame, to be
%                        placed into parameters.hzeeman and sent
%                        to the pulse sequence
%
%                        'iso_eq' - thermal equilibrium is com-
%                        computed using the isotropic part of 
%                        the Hamiltonian, and sent to the pulse
%                        sequence via parameters.rho0
%
%                        'aniso_eq' - thermal equilibrium is re-
%                        computed using the full anisotropic Ha-
%                        miltonian at each orientation, and sent
%                        to pulse sequence via parameters.rho0
%
%   parameters.rho0    - initial state; may be a function handle
%                        that depends on the three Euler angles
%                        in ZYZ active convention
%  
%   parameters.serial  - if set to true, disables automatic pa-
%                        rallelisation
%
%   parameters.sum_up  - if set to false, causes the pulse sequ-
%                        ence output at each orientation to be
%                        returned instead of the powder average
%
%   parameters.*       - additional subfields may be required by
%                        the pulse sequence - check its documen-
%                        tation page 
% 
%   assumptions        - context-specific assumptions ('nmr', 'epr',
%                        'labframe', etc.) - see the pulse sequence
%                        header for information on this setting.
%
% Outputs:
%
%   answer    - powder average of whatever it is that the pulse 
%               sequence returns; if parameters.sum_up is set to
%               false, a cell array of outputs at each orienta-
%               tion is returned
%
%   sph_grid -  powder averaging grid data structure with three 
%               Euler angles and weights for each point
%
% Note: THIS IS FOR STATIC POWDERS - use singlerot for MAS simulations.
%
% Note: arbitrary order rotating frame transformation is supported, inc-
%       luding infinite order. See the header of rotframe.m for further
%       information.
%
% Note: the function supports parallel processing via Matlab's Distri-
%       buted Computing Toolbox - different system orientations are eva-
%       luated on different labs.
%
% ledwards@cbs.mpg.de
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=powder.m>

function [answer,sph_grid]=powder(spin_system,pulse_sequence,...
                                  parameters,assumptions)

% Show the banner
banner(spin_system,'sequence_banner'); 

% Set common defaults
parameters=defaults(spin_system,parameters);

% Check consistency
grumble(spin_system,pulse_sequence,parameters,assumptions);

% Report to the user
report(spin_system,'building the Liouvillian...');

% Get the lab frame Zeeman operator if needed
if ismember('zeeman_op',parameters.needs)
    report(spin_system,'building the lab frame Zeeman operator...');
    [ZI,ZQ]=hamiltonian(assume(spin_system,'labframe','zeeman'));
else
    ZI=[]; ZQ=[]; % Work around the parfor bug
end

% Get the lab frame Hamiltonian if needed
if ismember('iso_eq',parameters.needs)
    
    % Isotropic part of the lab frame Hamiltonian
    report(spin_system,'building the lab frame Hamiltonian...');
    HL=hamiltonian(assume(spin_system,'labframe'),'left'); 

    % Compute thermal equilibrium here
    parameters.rho0=equilibrium(spin_system,HL); QL=[];

elseif ismember('aniso_eq',parameters.needs)
    
    % Isotropic + anisotropic part and thermal equilibrium later
    report(spin_system,'building the lab frame Hamiltonian...');
    [HL,QL]=hamiltonian(assume(spin_system,'labframe'),'left');

else

    % Work around the parfor bug
    HL=[]; QL=[]; 

end

% Set the assumptions
spin_system=assume(spin_system,assumptions);

% Get the Hamiltonian
[I,Q]=hamiltonian(spin_system);

% Get kinetics superoperator
K=kinetics(spin_system);

% Add offsets to the isotropic part
I=frqoffset(spin_system,I,parameters);

% Get carrier operators
C=cell(size(parameters.rframes));
for n=1:numel(parameters.rframes)
    C{n}=carrier(spin_system,parameters.rframes{n}{1});
end

% Get problem dimensions
parameters.spc_dim=1; parameters.spn_dim=size(I,1);

% Get the averaging grid as a structure
sph_grid=load([spin_system.sys.root_dir '/kernel/grids/' ...
               parameters.grid],'alphas','betas','gammas','weights');
           
% Assign local variables
alphas=sph_grid.alphas; betas=sph_grid.betas; 
gammas=sph_grid.gammas; weights=sph_grid.weights;
n_orients=numel(weights);

% Preallocate answer array
ans_array=cell(n_orients,1);

% Run serially if needed
if isfield(parameters,'serial')&&...
           parameters.serial
       
    % Serial execution
    nworkers=0;
    
    % Inform the user
    report(spin_system,'WARNING: parallelisation turned off by the user.');
    
else
    
    % Parallel execution
    nworkers=min([poolsize n_orients]);

end

% Inform the user and silence the output
prev_setting=spin_system.sys.output; ss_prev.sys.output=spin_system.sys.output;
report(spin_system,['powder simulation with ' num2str(n_orients) ' orientations.']);
if ~isfield(parameters,'verbose')||(parameters.verbose==0)
    report(spin_system,'pulse sequence silenced to avoid excessive output.')
    spin_system.sys.output='hush';
end

% Parfor rigging
if ~isworkernode
    D=parallel.pool.DataQueue;
    afterEach(D,@(~)parfor_progr);
    orients_done=0; last_toc=0;
    tic; ticBytes(gcp); do_diag=true;
else
    do_diag=false; D=[];
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
parfor (n=1:n_orients,nworkers)
    
    % Localise parameter array
    localpar=parameters;

    % Pass the current orientation to the pulse sequence
    localpar.current_angles=[alphas(n) betas(n) gammas(n)];
    
    % Get the full Hamiltonian at the current orientation
    H=I+orientation(Q,[alphas(n) betas(n) gammas(n)]); H=(H+H')/2;
    
    % Get the lab frame Zeeman operator at the current orientation
    if ismember('zeeman_op',parameters.needs)
        Z=ZI+orientation(ZQ,[alphas(n) betas(n) gammas(n)]);
        localpar.hzeeman=(Z+Z')/2;
    end
    
    % Compute initial state at the current orientation
    if ismember('aniso_eq',parameters.needs)

        % Anisotropic thermal equilibrium state from Hamiltonian and temperature
        localpar.rho0=equilibrium(spin_system,HL,QL,[alphas(n) betas(n) gammas(n)]);

    elseif isfield(parameters,'rho0')&&isa(parameters.rho0,'function_handle')

        % Anisotropic initial condition specified by the user
        rho_init=parameters.rho0; localpar.rho0=rho_init(alphas(n),betas(n),gammas(n));

    end
    
    % Apply rotating frames
    for k=1:numel(localpar.rframes)
        H=rotframe(spin_system,C{k},H,localpar.rframes{k}{1},...
                                      localpar.rframes{k}{2}); %#ok<PFBNS>
    end
    
    % Get the relaxation superoperator at the current orientation
    R=relaxation(spin_system,[alphas(n) betas(n) gammas(n)]);
    
    % Report to the user
    report(spin_system,'running the pulse sequence...');
    
    % Run the simulation (it might return a structure)
    ans_array{n}=pulse_sequence(spin_system,localpar,H,R,K); %#ok<PFBNS>

    % Report progress
    if do_diag, send(D,n); end
    
end

% Decide the return array
if parameters.sum_up
    
    % Return weighted sum
    answer=sph_grid.weights(1)*ans_array{1};
    for n=2:numel(ans_array)
        answer=answer+sph_grid.weights(n)*ans_array{n};
    end
    
    % Inform the user
    report(spin_system,'returning powder averaged pulse sequence output...');
    
else
    
    % Return components
    answer=ans_array;
    
    % Inform the user
    report(spin_system,'returning pulse sequence outputs at each orientation...');
    
end

% Unsilence the output
spin_system.sys.output=prev_setting;

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
    report(spin_system,'parameters.rframes field not set, assuming no additional rotating frames.');
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
    report(spin_system,'parameters.sum_up field not set, returning grid integral.');
    parameters.sum_up=1;
end
if ~isfield(parameters,'needs')
    report(spin_system,'parameters.needs field not set, assumpting empty.');
    parameters.needs={};
end
end

% Consistency checking
function grumble(spin_system,pulse_sequence,parameters,assumptions)

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

% Additional needs
if isfield(parameters,'needs')
    if ismember('iso_eq',parameters.needs)&&ismember('aniso_eq',parameters.needs)
        error('iso_eq and aniso_eq needs cannot be specified simultaneously.');
    end
    if isfield(parameters,'rho0')
        if ismember('iso_eq',parameters.needs)||...
           ismember('aniso_eq',parameters.needs)
            error('parameters.needs cannot request initial condition when parameters.rho0 is specified.');
        end
    end
end

end

% Do not confuse altruism with kindness. The irreducible primary of
% altruism is self-sacrifice - which means self-immolation, self-
% abnegation, self-denial. Do not hide behind such superficialities
% as whether you should or should not give a dime to a beggar. This
% is not the issue. The issue is whether you do or do not have the
% right to exist without giving him that dime. The issue is whether
% the need of others is the first mortgage on your life and the mo-
% ral purpose of your existence. Any man of self-esteem will answer:
% No. Altruism says: Yes.
%
% Ayn Rand

