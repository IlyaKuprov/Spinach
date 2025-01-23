% Magic angle spinning DNP simulation, returning the rotor period 
% averaged steady state magnetization. This function takes a lot 
% of inspiration from the code donated by Fred Mentink, please ci-
% te Fred's papers if you are using it. Syntax:
%
%                dnp=masdnp(spin_system,parameters)
%
% Parameters:
%
%     parameters.spins     -  the spins to microwave
%
%     parameters.rate      -  spinning rate, Hz
% 
%     parameters.axis      -  spinning axis direction vector.
%
%     parameters.max_rank  -  rotor discretization grid rank
%
%     parameters.mw_pwr    -  microwave power, rad/s
%
%     parameters.mw_frq    -  microwave frequency, Hz
%
%     parameters.mw_time   -  microwave irradiation duration 
%                             before the average magnetistion
%                             is computed, seconds
%
%     parameters.grid      -  the name of the spherical avera-
%                             ging grid
%
%     parameters.coil      -  detection state
%
%     parameters.verbose   -  set this to 1 to enable diag-
%                             nostic output
%
% Outputs:
%
%     dnp - enhancement of the user-specified state relative to
%           the thermal equilibrium
%
% Note: increase the rotor rank and the spherical grid size until
%       the answer stops changing. You will likely need huge values
%       for both parameters.
%
% Note: this function must be called directly, without a context
%       wrapper.
%
% ilya.kuprov@weizmann.ac.il
% fmentink@magnet.fsu.edu
%
% <https://spindynamics.org/wiki/index.php?title=masdnp.m>

function dnp=masdnp(spin_system,parameters)

% Check consistency
grumble(spin_system,parameters);

% Convert microwave frequency into offset
parameters.offset=parameters.mw_frq-spin_system.inter.magnet*spin(parameters.spins{1})/(2*pi);

% Microwave operator
Hmw=(operator(spin_system,'L+',parameters.spins{1})+...
     operator(spin_system,'L-',parameters.spins{1}))/2;
 
% Relaxation superoperator
R=relaxation(spin_system);

% Thermal equilibrium state
rho_eq=equilibrium(spin_system,hamiltonian(assume(spin_system,'labframe'),'left'));

% Get the averaging grid
sph_grid=load([spin_system.sys.root_dir '/kernel/grids/' parameters.grid '.mat']);

% Shut up and inform the user
report(spin_system,['powder average being computed over ' ...
                    num2str(numel(sph_grid.weights)) ' orientations.']);
if ~parameters.verbose
    report(spin_system,'pulse sequence silenced to avoid excessive output.');
    spin_system.sys.output='hush';
end

% Get the answer going
dnp=0;

% Loop over the grid weights
parfor n=1:numel(sph_grid.weights) %#ok<*PFBNS>
    
    % Set the current orientation
    localpar=parameters;
    localpar.orientation=[sph_grid.alphas(n) ...
                          sph_grid.betas(n)  ...
                          sph_grid.gammas(n)];      
    localpar.masframe='rotor'; localpar.rframes={};
    
    % Rotor stack generation
    H=rotor_stack(spin_system,localpar,'esr');
    nsteps=numel(H); dt=1/(nsteps*localpar.rate);
    
    % Rotor period integration
    nsteps=numel(H); P=speye(size(H{1}));
    for k=1:nsteps
        P=propagator(spin_system,H{k}+localpar.mw_pwr*Hmw+1i*R,dt)*P; 
    end
    
    % Effective rotor period Liouvillian
    L_eff=1i*localpar.rate*logm(P);

    % Evolve for the equilibration time
    rho_st=evolution(spin_system,L_eff,[],rho_eq,parameters.mw_time,1,'final');

    % Rotor period trajectory
    rho=zeros([numel(rho_eq) nsteps],'like',1i); rho(:,1)=rho_st;
    for k=2:nsteps
        rho(:,k)=step(spin_system,H{k}+localpar.mw_pwr*Hmw+1i*R,rho(:,k-1),dt);
    end

    % Enhancement factor
    Hz_dnp=mean(localpar.coil'*rho);
    Hz_eq=localpar.coil'*rho_eq;
    dnp=dnp+sph_grid.weights(n)*Hz_dnp/Hz_eq;
    
end
        
end

% Consistency enforcement
function grumble(spin_system,parameters)
if ~ismember(spin_system.bas.formalism,{'zeeman-liouv','sphten-liouv'})
    error('this function is only available in Liouville space.');
end
if ~isfield(parameters,'max_rank')
    error('parameters.max_rank subfield must be present.');
elseif (~isnumeric(parameters.max_rank))||(~isreal(parameters.max_rank))||...
       (mod(parameters.max_rank,1)~=0)||(parameters.max_rank<0)
    error('parameters.max_rank must be a positive real integer.');
end
if ~isfield(parameters,'rate')
    error('parameters.rate subfield must be present.');
elseif (~isnumeric(parameters.rate))||(~isreal(parameters.rate))
    error('parameters.rate must be a real number.');
end
if ~isfield(parameters,'axis')
    error('parameters.axis subfield must be present.');
elseif (~isnumeric(parameters.axis))||(~isreal(parameters.axis))||...
       (~isrow(parameters.axis))||(numel(parameters.axis)~=3)
    error('parameters.axis must be a row vector of three real numbers.');
end
if isempty(parameters.spins)
    error('parameters.spins variable cannot be empty.');
elseif ~iscell(parameters.spins)
    error('parameters.spins variable must be a cell array.');
elseif ~all(cellfun(@ischar,parameters.spins))
    error('all elements of parameters.spins cell array must be strings.');
elseif any(~ismember(parameters.spins,spin_system.comp.isotopes))
    error('parameters.spins refers to a spin that is not present in the system.');
end
if ~isfield(parameters,'grid')
    error('spherical averaging grid must be specified in parameters.grid variable.');
elseif isempty(parameters.grid)
    error('parameters.grid variable cannot be empty.');
elseif ~ischar(parameters.grid)
    error('parameters.grid variable must be a character string.');
end
if ~isfield(parameters,'coil')
    error('detection state must be specified in parameters.coil variable.');
end
if ~isfield(parameters,'mw_pwr')
    error('microwave power must be specified in parameters.mw_pwr field.');
end
if (~isnumeric(parameters.mw_pwr))||(~isreal(parameters.mw_pwr))||...
   (numel(parameters.mw_pwr)~=1)||(parameters.mw_pwr<0)
    error('parameters.mw_pwr must be a non-negative real scalar.');
end
if ~isfield(parameters,'mw_frq')
    error('microwave frequency must be specified in parameters.mw_frq field.');
end
if (~isnumeric(parameters.mw_frq))||(~isreal(parameters.mw_frq))||...
   (numel(parameters.mw_frq)~=1)
    error('parameters.mw_frq must be a real scalar.');
end
if ~isfield(parameters,'mw_time')
    error('microwaving time must be specified in parameters.mw_time field.');
end
if (~isnumeric(parameters.mw_time))||(~isreal(parameters.mw_time))||...
   (numel(parameters.mw_time)~=1)||(parameters.mw_time<0)||...
   (mod(parameters.mw_time,1)~=0)
    error('parameters.mw_time must be a non-negative real integer.');
end
end

% If you have a garden and a library, you have 
% everything you need. 
%
% Marcus Tullius Cicero

