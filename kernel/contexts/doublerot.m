% Fokker-Planck double angle spinning context. Generates a Liouvillian
% superoperator and passes it on to the pulse sequence function, which
% should be supplied as a handle. Syntax:
%
%      [answer,sph_grid]=doublerot(spin_system,pulse_sequence,...
%                                  parameters,assumptions)
%
% where pulse sequence is a function handle to one of the pulse sequences
% located in the experiments directory, assumptions is a string that would
% be passed to assume.m when the Hamiltonian is built and parameters is a
% structure with the following subfields:
%
%   parameters.rate_outer - outer rotor spinning rate in Hz
%
%   parameters.rate_inner - inner rotor spinning rate in Hz
%
%   parameters.axis_outer - spinning axis of the outer rotor,
%                           given as a normalized 3-element
%                           vector
%
%   parameters.axis_inner - spinning axis of the inner rotor,
%                           given as a normalized 3-element
%                           vector
%
%   parameters.rank_outer - maximum harmonic rank to retain in
%                           the solution (increase till conver-
%                           gence is achieved, approximately
%                           equal to the number of spinning si-
%                           debands in the spectrum) for the
%                           outer rotor
%
%   parameters.rank_inner - maximum harmonic rank to retain in
%                           the solution (increase till conver-
%                           gence is achieved, approximately
%                           equal to the number of spinning si-
%                           debands in the spectrum) for the
%                           inner rotor
%
%   parameters.rframes  - rotating frame specification, e.g.
%                         {{'13C',2},{'14N,3}} requests second
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
%   parameters.serial   - if set to true, disables automatic pa-
%                         rallelisation
%
%   parameters.sum_up   - when set to 1 (default), returns the
%                         powder average. When set to 0, returns
%                         individual answers for each point in 
%                         the powder as a cell array.
%
%   parameters.*        - additional subfields may be required by your
%                         pulse sequence - check its documentation page 
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
% Note: the state projector assumes a powder -- single crystal DOR is not
%       currently supported.
%
% Note: the function supports parallel processing via Matlab's Distri-
%       buted Computing Toolbox - different system orientations are eva-
%       luated on different labs.
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=doublerot.m>

function [answer,sph_grid]=doublerot(spin_system,pulse_sequence,...
                                     parameters,assumptions)

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

% Count rotor trajectory points
npoints_inner=2*parameters.rank_inner+1;
npoints_outer=2*parameters.rank_outer+1;
npoints_total=npoints_outer*npoints_inner;

% Get problem dimensions
spc_dim=npoints_total; spn_dim=size(H,1);
report(spin_system,['lab space problem dimension     ' num2str(spc_dim)]);
report(spin_system,['spin space problem dimension    ' num2str(spn_dim)]);
report(spin_system,['Fokker-Planck problem dimension ' num2str(spc_dim*spn_dim)]);
parameters.spc_dim=spc_dim; parameters.spn_dim=spn_dim;

% Compute spectral derivative operators
[traj_inner,d_dphi_inner]=fourdif(npoints_inner,1);
[traj_outer,d_dphi_outer]=fourdif(npoints_outer,1);

% Compute rotor phase tracks
phases_outer=kron(traj_outer,ones(npoints_inner,1));
phases_inner=kron(ones(npoints_outer,1),traj_inner);

% Compute double spinning operator
M=2*pi*kron((parameters.rate_outer*kron(d_dphi_outer,speye(size(d_dphi_inner)))+...
             parameters.rate_inner*kron(speye(size(d_dphi_outer)),d_dphi_inner)),speye(size(H)));

% Get rotor axis orientations
[phi_inner,theta_inner,~]=cart2sph(parameters.axis_inner(1),...
                                   parameters.axis_inner(2),...
                                   parameters.axis_inner(3));
[phi_outer,theta_outer,~]=cart2sph(parameters.axis_outer(1),...
                                   parameters.axis_outer(2),...
                                   parameters.axis_outer(3));
theta_inner=pi/2-theta_inner; theta_outer=pi/2-theta_outer;

% Get carrier operators
C=cell(size(parameters.rframes));
for n=1:numel(parameters.rframes)
    C{n}=carrier(spin_system,parameters.rframes{n}{1});
end

% Get relaxation and kinetics 
R=relaxation(spin_system); K=kinetics(spin_system);

% Get the averaging grid
sph_grid=load([spin_system.sys.root_dir '/kernel/grids/' parameters.grid '.mat']);

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

% Powder averaged spectrum
parfor (q=1:numel(weights),nworkers) %#ok<*PFBNS>
    
    % Preallocate Liouvillian blocks
    L=cell(npoints_total,npoints_total);
    for n=1:npoints_total
        for k=1:npoints_total
            L{n,k}=krondelta(n,k)*H;
        end
    end
    
    % Build Liouvillian blocks
    for n=1:npoints_total
        
        % Loop over spherical ranks
        for r=1:numel(Q)
            
            % Compute outer rotor axis tilt
            D_lab2out=wigner(r,phi_outer,theta_outer,0);

            % Compute inner rotor axis tilt
            D_out2inn=wigner(r,phi_inner,theta_inner,0);
    
            % Compute crystallite orientation
            D_mol2inn=wigner(r,alphas(q),betas(q),gammas(q));
            
            % Compute outer rotor rotation
            D_outer=wigner(r,0,0,phases_outer(n));
            
            % Compute inner rotor rotation
            D_inner=wigner(r,0,0,phases_inner(n));
            
            % Compose rotations
            D=D_lab2out*D_outer*D_out2inn*D_inner*D_mol2inn;
            
            % Build the block
            for k=1:(2*r+1)
                for m=1:(2*r+1)
                    L{n,n}=L{n,n}+D(k,m)*Q{r}{k,m};
                end
            end
            
        end
        
        % Apply rotating frames
        for k=1:numel(parameters.rframes)
            L{n,n}=rotframe(spin_system,C{k},(L{n,n}+L{n,n}')/2,parameters.rframes{k}{1},parameters.rframes{k}{2});
        end
        
    end
    
    % Assemble the liouvillian
    L=clean_up(spin_system,cell2mat(L),spin_system.tols.liouv_zero);
    
    % Report to the user
    report(spin_system,'running the pulse sequence...');
    
    % Run the pulse sequence
    ans_array{q}=pulse_sequence(spin_system,parameters,L+1i*M,R,K);
    
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
    
    % Return components
    answer=ans_array;
   
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

% Rotor ranks
if ~isfield(parameters,'rank_outer')
    error('parameters.rank_outer subfield must be present.');
elseif (~isnumeric(parameters.rank_outer))||(~isreal(parameters.rank_outer))||...
       (mod(parameters.rank_outer,1)~=0)||(parameters.rank_outer<0)
    error('parameters.rank_outer must be a positive real integer.');
end
if ~isfield(parameters,'rank_inner')
    error('parameters.rank_inner subfield must be present.');
elseif (~isnumeric(parameters.rank_inner))||(~isreal(parameters.rank_inner))||...
       (mod(parameters.rank_inner,1)~=0)||(parameters.rank_inner<0)
    error('parameters.rank_inner must be a positive real integer.');
end

% Spinning rates
if ~isfield(parameters,'rate_outer')
    error('parameters.rate_outer subfield must be present.');
elseif (~isnumeric(parameters.rate_outer))||(~isreal(parameters.rate_outer))
    error('parameters.rate_outer must be a real number.');
end
if ~isfield(parameters,'rate_inner')
    error('parameters.rate_inner subfield must be present.');
elseif (~isnumeric(parameters.rate_inner))||(~isreal(parameters.rate_inner))
    error('parameters.rate_inner must be a real number.');
end

% Spinning axes
if ~isfield(parameters,'axis_outer')
    error('parameters.axis_outer subfield must be present.');
elseif (~isnumeric(parameters.axis_outer))||(~isreal(parameters.axis_outer))||...
       (~isrow(parameters.axis_outer))||(numel(parameters.axis_outer)~=3)
    error('parameters.axis_outer must be a row vector of three real numbers.');
end
if ~isfield(parameters,'axis_inner')
    error('parameters.axis_inner subfield must be present.');
elseif (~isnumeric(parameters.axis_inner))||(~isreal(parameters.axis_inner))||...
       (~isrow(parameters.axis_inner))||(numel(parameters.axis_inner)~=3)
    error('parameters.axis_inner must be a row vector of three real numbers.');
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

% Show me a hero and I'll write you a tragedy.
%
% F. Scott Fitzgerald

