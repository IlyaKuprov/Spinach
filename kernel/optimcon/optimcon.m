% Validates optimal control options and updates the spin system
% object. Syntax:
% 
%          spin_system=optimcon(spin_system,control)
%
% Parameters:
%
%     spin_system  - primary Spinach data structure,
%                    created by create.m and updated
%                    by basis.m functions
%
%     control      - control data structure described 
%                    in detail in the online manual
%
% Outputs:
%
%     spin_system  - updated Spinach data structure
%
% david.goodwin@inano.au.dk
% u.rasulov@soton.ac.uk
% ilya.kuprov@weizmann.ac.il
% m.keitel@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=optimcon.m>

function spin_system=optimcon(spin_system,control)

% Consistency check
grumble(spin_system,control)

% Delete the previous control structure
if isfield(spin_system,'control')
    spin_system=rmfield(spin_system,'control');
    report(spin_system,'WARNING: pre-existing control structure deleted.');
end

% Show the banner 
banner(spin_system,'optimcon');

% Process fidelity type
if isfield(control,'fidelity')

    % Input validation
    if ~ischar(control.fidelity)
        error('control.fidelity must be a character string.');
    end
    
    % Absorb the input
    spin_system.control.fidelity=control.fidelity;
    control=rmfield(control,'fidelity');

else

    % Default is Re(<targ|P|init>)
    spin_system.control.fidelity='real';

end

% Inform the user
switch spin_system.control.fidelity
    
    case 'real'
        
        % Real part of the overlap
        report(spin_system,[pad('Fidelity measure, range [-1,+1]',60)...
                            pad('Re(<target|rho(T)>)',20)]);
                          
    case 'imag'
        
        % Imaginary part of the overlap
        report(spin_system,[pad('Fidelity measure, range [-1,+1]',60)...
                            pad('Im(<target|rho(T)>)',20)]);
                          
    case 'square'
        
        % Absolute square of the overlap
        report(spin_system,[pad('Fidelity measure, range [0,+1]',60)...
                            pad('|<target|rho(T)>|^2',20)]);
                        
    otherwise
        
        % Complain and bomb out
        error('control.fidelity can be ''real'', ''imag'', or ''square''.');
                          
end

% Process integrator type
if isfield(control,'integrator')

    % Input validation
    if ~ischar(control.integrator)
        error('control.integrator must be a character string.');
    end
    if ~ismember(control.integrator,{'rectangle','trapezium'})
        error('control.intgrator must be either ''rectangle'' or ''trapezium''.');
    end

    % Absorb integrator type
    spin_system.control.integrator=control.integrator;
    control=rmfield(control,'integrator');

else

    % Default is piecewise-constant
    spin_system.control.integrator='rectangle';

end

% Inform the user
report(spin_system,[pad('Equation of motion integrator',60) ...
                    spin_system.control.integrator]);

% Process control operators
if isfield(control,'operators')

    % Input validation
    if (~iscell(control.operators))||(~all(cellfun(@ismatrix,control.operators(:))))
        error('control.operators must be a cell array of matrices.');
    end

    % In Hilbert space, disallow non-Hermitian controls
    if strcmp('zeeman-hilb',spin_system.bas.formalism)

        % Over controls
        for n=1:numel(control.operators)

            % Get pertinent norms
            norm_a=cheap_norm(control.operators{n}-...
                              control.operators{n}');
            norm_b=cheap_norm(control.operators{n});

            % Check the norms
            if norm_a>1e-10*norm_b
                error('all control generators must be Hermitian in Hilbert space.');
            end

        end

    end

    % Control operator count
    spin_system.control.ncontrols=numel(control.operators);

    % Absorb control operators and clean them up
    for n=1:spin_system.control.ncontrols
        spin_system.control.operators{n}=clean_up(spin_system,control.operators{n},...
                                                  spin_system.tols.liouv_zero);
    end
    control=rmfield(control,'operators');

    % Inform the user
    report(spin_system,[pad('Number of control operators',60) ...
                        int2str(spin_system.control.ncontrols)]);

    % Precompute control-control commutators if necessary
    if strcmp(spin_system.control.integrator,'trapezium')

        % Commutator array and index preallocation
        spin_system.control.cc_comm=cell(spin_system.control.ncontrols, ...
                                         spin_system.control.ncontrols);
        spin_system.control.cc_comm_idx=false(spin_system.control.ncontrols, ...
                                              spin_system.control.ncontrols);

        % Commutator calculation
        for n=1:spin_system.control.ncontrols
            for m=1:spin_system.control.ncontrols

                % Commutator itself (usually very sparse, not a pain to store)
                spin_system.control.cc_comm{n,m}=comm(spin_system.control.operators{n},...
                                                      spin_system.control.operators{m});

                % Keep a logical table of commutators (true when they commute)
                comm_norm=norm(spin_system.control.cc_comm{n,m},1);
                spin_system.control.cc_comm_idx(n,m)=(comm_norm<spin_system.tols.liouv_zero);

            end
        end

        % Inform the user
        n_commute=(nnz(spin_system.control.cc_comm_idx)-numel(spin_system.control.operators))/2;
        report(spin_system,[pad('Number of commuting control pairs',60) int2str(n_commute)]);

    end
    
else

    % Complain and bomb out
    error('control operators must be supplied in control.operators field.');

end

% Process initial states
if isfield(control,'rho_init')
    
    % Input validation
    if ~iscell(control.rho_init)
        error('control.rho_init must be a cell array of vectors.');
    end

    % In-depth checking
    for n=1:numel(control.rho_init)

        switch spin_system.bas.formalism

            case {'sphten-liouv','zeeman-liouv'}

                if (~isnumeric(control.rho_init{n}))||(~iscolumn(control.rho_init{n}))
                    error('control.rho_init must be a cell array of column vectors.');
                end

            case 'zeeman-hilb'

                if (~isnumeric(control.rho_init{n}))||(size(control.rho_init{n},1)~=...
                                                       size(control.rho_init{n},2))
                    error('control.rho_init must be a cell array of square matrices.');
                end

            case 'zeeman-wavef'

                if (~isnumeric(control.rho_init{n}))||(size(control.rho_init{n},2)~=1)
                    error('control.rho_init must be a cell array of column vectors.');
                end

            otherwise

                error('unrecognised formalism specification.');

        end

    end
    
    % Absorb initial states
    spin_system.control.rho_init=control.rho_init; 
    control=rmfield(control,'rho_init');

    % Inform the user
    report(spin_system,[pad('Initial states per ensemble member',60) ...
                        int2str(numel(spin_system.control.rho_init))]);

else

    % Complain and bomb out
    error('initial states must be supplied in control.rho_init field.');

end

% Process target states
if isfield(control,'rho_targ')

    % Input validation
    if ~iscell(control.rho_targ)
        error('control.rho_targ must be a cell array of vectors.');
    end
    
    % In-depth checking
    for n=1:numel(control.rho_targ)

        switch spin_system.bas.formalism

            case {'sphten-liouv','zeeman-liouv'}

                if (~isnumeric(control.rho_targ{n}))||(~iscolumn(control.rho_targ{n}))
                    error('control.rho_targ must be a cell array of column vectors.');
                end

            case 'zeeman-hilb'

                if (~isnumeric(control.rho_targ{n}))||(size(control.rho_targ{n},1)~=...
                                                       size(control.rho_targ{n},2))
                    error('control.rho_targ must be a cell array of square matrices.');
                end

            case 'zeeman-wavef'

                if (~isnumeric(control.rho_targ{n}))||(size(control.rho_targ{n},2)~=1)
                    error('control.rho_targ must be a cell array of column vectors.');
                end

            otherwise

                error('unrecognised formalism specification.');

        end

    end

    % Source-target pair consistency
    if numel(control.rho_targ)~=numel(spin_system.control.rho_init)
        error('control.rho_targ must have the same size as control.rho_init');
    end

    % Absorb target states
    spin_system.control.rho_targ=control.rho_targ; 
    control=rmfield(control,'rho_targ');

    % Inform the user
    report(spin_system,[pad('Target states per ensemble member',60) ...
                        int2str(numel(spin_system.control.rho_targ))]);

else

    % Complain and bomb out
    error('target states must be supplied in control.rho_targ field.');

end

% Stroboscopic steady state switch
if isfield(control,'steady')

    % Input validation
    if ~islogical(control.steady) || ~isscalar(control.steady)
        error('control.steady must be true() or false()');
    end

    % Absorb stroboscopic steady state switch
    spin_system.control.steady=control.steady; 
    control=rmfield(control,'steady');

    % Inform the user
    if spin_system.control.steady
        report(spin_system,[pad('Stroboscopic steady state',60)  'on']);
    else
        report(spin_system,[pad('Stroboscopic steady state',60)  'off']);
    end

else

    % Default is user-specified states
    spin_system.control.steady=false();

    % Inform the user
    report(spin_system,[pad('Stroboscopic steady state',60) 'off']);

end

% Process prefix function
if isfield(control,'prefix')

    % Input validation
    if ~isa(control.prefix,'function_handle')
        error('control.prefix must contain a function handle.');
    end
    
    % Absorb prefix function
    spin_system.control.prefix=control.prefix;
    control=rmfield(control,'prefix');
    
    % Inform the user
    report(spin_system,[pad('Prefix sequence function',60) '@' ...
                        char(spin_system.control.prefix)]);
                    
else
    
    % Empty prefix function
    spin_system.control.prefix=[];
    
    % Inform the user
    report(spin_system,[pad('Prefix sequence function',60) 'none']);
    
end

% Process suffix function
if isfield(control,'suffix')

    % Input validation
    if ~isa(control.suffix,'function_handle')
        error('control.suffix must contain a function handle.');
    end
    
    % Absorb suffix function
    spin_system.control.suffix=control.suffix;
    control=rmfield(control,'suffix');
    
    % Inform the user
    report(spin_system,[pad('Suffix sequence function',60) '@' ...
                        char(spin_system.control.suffix)]);
                    
else
    
    % Empty suffix function
    spin_system.control.suffix=[];
    
    % Inform the user
    report(spin_system,[pad('Suffix sequence function',60) 'none']);
    
end

% Process distortion functions
if isfield(control,'distortion')

    % Disallow Newton-Raphson
    if ismember(control.method,{'newton','goodwin'})
        error('distortion handling not implemented for Newton-Raphson methods.');
    end

    % Input validation
    if ~iscell(control.distortion)
        error('control.distortion must be a cell array of function handles.');
    end
    for n=1:numel(control.distortion)
        if ~isa(control.distortion{n},'function_handle')
            error('control.distortion must be a cell array of function handles.');
        end
    end
    
    % Absorb distortion functions
    spin_system.control.distortion=control.distortion;
    control=rmfield(control,'distortion');
    
    % Inform the user
    report(spin_system,[pad('Waveform distortion stage count',60)   ...
                        int2str(size(spin_system.control.distortion,2))]);
    report(spin_system,[pad('Waveform distortion ensemble size',60) ...
                        int2str(size(spin_system.control.distortion,1))]);

    % Block methods that use exact Hessians
    if ismember(control.method,{'newton','goodwin'})
        error('waveform distortions are only available with LBFGS optimiser.');
    end
                    
else
    
    % Only one function that does nothing
    spin_system.control.distortion={@no_dist};
    
    % Inform the user
    report(spin_system,[pad('Waveform distortion stages',60) 'none']);
    report(spin_system,[pad('Waveform distortion ensemble',60) 'none']);
    
end

% Process dead time
if isfield(control,'dead_time')

    % Input validation
    if (~isnumeric(control.dead_time))||(~isreal(control.dead_time))||...
       (~isscalar(control.dead_time))||(control.dead_time<0)
        error('control.dead_time must be a non-negative real number.');
    end
    
    % Absorb dead time
    spin_system.control.dead_time=control.dead_time;
    control=rmfield(control,'dead_time');
    
else
    
    % Default is zero dead time
    spin_system.control.dead_time=0;
    
end

% Inform the user
report(spin_system,[pad('Dead time, seconds',60) ...
                         num2str(spin_system.control.dead_time,'%.9g')]);

% Process power levels
if isfield(control,'pwr_levels')

    % Input validation
    if (~isnumeric(control.pwr_levels))||(~isreal(control.pwr_levels))||...
       (~isrow(control.pwr_levels))||any(control.pwr_levels(:)<=0)
        error('control.pwr_levels must be a row vector of positive real numbers.');
    end

    % Absorb power levels
    spin_system.control.pwr_levels=control.pwr_levels; 
    control=rmfield(control,'pwr_levels');

    % Inform the user
    report(spin_system,[pad('Control power ensemble size',60) ...
                        int2str(numel(spin_system.control.pwr_levels))]);
    report(spin_system,[pad('   Max power multiplier, Hz',60) ...
                        num2str(max(spin_system.control.pwr_levels)/(2*pi),'%.9g')]);
    report(spin_system,[pad('   Min power multiplier, Hz',60) ...
                        num2str(min(spin_system.control.pwr_levels)/(2*pi),'%.9g')]);

else

    % Complain and bomb out
    error('power levels must be specified in control.pwr_levels field.');

end

% Process offset distributions
if isfield(control,'offsets')&&isfield(control,'off_ops')
    
    % Basic type validation
    if ~iscell(control.offsets)
        error('control.offsets must be a cell array of row vectors.');
    end
    if ~iscell(control.off_ops)
        error('control.offsets must be a cell array of square matrices.');
    end
    if numel(control.offsets)~=numel(control.off_ops)
        error('control.offsets and control.off_ops must have the same size.');
    end

    % Content validation
    for n=1:numel(control.offsets)
        if (~isnumeric(control.offsets{n}))||...
           (~isreal(control.offsets{n}))||...
           (~isrow(control.offsets{n}))
            error('elements of control.offsets must be row vectors of real numbers.');
        end
        if (~isnumeric(control.off_ops{n}))||...
           (size(control.off_ops{n},1)~=size(control.off_ops{n},2))
            error('elements of control.off_ops must be square matrices.');
        end
    end

    % Absorb offset distributions
    spin_system.control.offsets=control.offsets; 
    control=rmfield(control,'offsets');

    % Clean up and absorb offset operators
    spin_system.control.off_ops=clean_up(spin_system,control.off_ops,...
                                         spin_system.tols.liouv_zero);
    control=rmfield(control,'off_ops');

    % Inform the user
    report(spin_system,[pad('Number of offset channels',60) ...
                        int2str(numel(spin_system.control.offsets))]);
    for n=1:numel(spin_system.control.offsets)
        report(spin_system,[pad(['   Channel ' int2str(n) ', ensemble size'],60) ...
                           int2str(numel(spin_system.control.offsets{n}))]);
        report(spin_system,[pad(['   Channel ' int2str(n) ', min offset, Hz'],60) ...
                           num2str(min(spin_system.control.offsets{n}),'%+.9g')]);
        report(spin_system,[pad(['   Channel ' int2str(n) ', max offset, Hz'],60) ...
                           num2str(max(spin_system.control.offsets{n}),'%+.9g')]);
    end
                
else
    
    % Default is no offsets
    spin_system.control.offsets={};
    spin_system.control.off_ops={};
    
end

% Process sequence timing
if isfield(control,'pulse_dt')

    % Input validation
    if (~isnumeric(control.pulse_dt))||(~isreal(control.pulse_dt))||...
       (~isrow(control.pulse_dt))||any(control.pulse_dt(:)<=0)
        error('control.pulse_dt must be a row vector of positive real numbers.');
    end

    % Extract timing parameters
    spin_system.control.pulse_dt=control.pulse_dt;
    spin_system.control.pulse_dur=sum(control.pulse_dt);
    spin_system.control.pulse_nsteps=numel(control.pulse_dt);
    control=rmfield(control,'pulse_dt');

    % Inform the user
    report(spin_system,[pad('Duration of the control sequence, seconds',60) ...
                        num2str(spin_system.control.pulse_dur,'%.9g')]);
    report(spin_system,[pad('Number of intervals in the control sequence',60)...
                        int2str(spin_system.control.pulse_nsteps)]);

    % Specify the number of values
    switch spin_system.control.integrator
        
        % Piecewise-constant
        case 'rectangle'
            
            % Each step is defined by one generator
            spin_system.control.pulse_ntpts=spin_system.control.pulse_nsteps;
        
        % Piecewise-linear
        case 'trapezium'

            % Each step is defined by left and right generators
            spin_system.control.pulse_ntpts=spin_system.control.pulse_nsteps+1;

        otherwise

            % Complain and bomb out
            error('unknown integrator type.');

    end

    % Inform the user
    report(spin_system,[pad('Control coefficient array length',60)...
                        int2str(spin_system.control.pulse_ntpts)]);

else

    % Complain and bomb out
    error('interval duration infomation in control.pulse_dt must be provided.');
  
end

% Process drift generators
if isfield(control,'drifts')

    % Input validation
    if ~iscell(control.drifts)
        error('control.drifts must be a cell array (over the ensemble) of cell arrays (over time) of matrices.');
    end
    for n=1:numel(control.drifts)
        if ~iscell(control.drifts{n})
            error('control.drifts must be a cell array (over ensemble) of cell arrays (over time) of matrices.');
        end
        if ~all(cellfun(@ismatrix,control.drifts{n}(:)))
            error('control.drifts must be a cell array (over ensemble) of cell arrays (over time) of matrices.');
        end
        if (numel(control.drifts{n})~=1)&&(numel(control.drifts{n})~=spin_system.control.pulse_ntpts)
            error(['need either 1 or ' int2str(spin_system.control.pulse_ntpts) ' elements in control.drift{' ...
                   int2str(n) '}, found ' int2str(numel(control.drifts{n})) ' elements.']);
        end
    end
    if ~all(cellfun(@numel,control.drifts(:))==numel(control.drifts{1}))
        error('all time-dependent drifts must have the same number of time slices.');
    end

    % In Hilbert space, disallow non-Hermitian drifts
    if strcmp('zeeman-hilb',spin_system.bas.formalism)
        
        % Over drift ensemble
        for n=1:numel(control.drifts)

            % Over time slices
            for k=1:numel(control.drifts{n})

                % Get pertinent norms
                norm_a=cheap_norm(control.drifts{n}{k}-...
                                  control.drifts{n}{k}');
                norm_b=cheap_norm(control.drifts{n}{k});

                % Check the norms
                if norm_a>1e-10*norm_b
                    error('all drift generators must be Hermitian in Hilbert space.');
                end

            end

        end

    end
    
    % Inform the user
    report(spin_system,[pad('Drift generator ensemble size',60) ...
                        int2str(numel(control.drifts))]);
    if isscalar(control.drifts{1})
        report(spin_system,[pad('Time-dependent drift generator',60) 'no']);
    else
        report(spin_system,[pad('Time-dependent drift generator',60) 'yes']);
    end

    % Put drift generators on pool ValueStore
    vs_labels=cell([numel(control.drifts) 1]);
    for n=1:numel(control.drifts)
        vs_labels{n}=['oc_drift_' num2str(n)];
    end
    store=gcp('nocreate').ValueStore;
    put(store,vs_labels,control.drifts);

    % Only keep the overall drift count
    spin_system.control.ndrifts=numel(control.drifts); 
    control=rmfield(control,'drifts');

else
    
    % Complain and bomb out
    error('drift generators must be supplied in control.drifts field.');

end

% Absorb the amplitude profile
if isfield(control,'amplitudes')

    % Store the profile (validated in grape_phase)
    spin_system.control.amplitudes=control.amplitudes;
    control=rmfield(control,'amplitudes');
 
end

% Process freeze mask
if isfield(control,'freeze')

    % Store the mask (validated in fmaxnewton)
    spin_system.control.freeze=logical(control.freeze);
    control=rmfield(control,'freeze');
    
    % Inform the user
    report(spin_system,[pad('Freeze mask supplied',60) 'yes']);
    
else

    % Store an empty array
    spin_system.control.freeze=[];
    
    % Inform the user
    report(spin_system,[pad('Freeze mask supplied',60) 'no']);
    
end      

% Process optimisation method
if isfield(control,'method')

    % Input validation
    if (~ischar(control.method))||(~ismember(control.method,{'lbfgs','newton','goodwin'}))
        error('control.method must be ''lbfgs'', ''newton'', or ''goodwin''.');
    end
    
    % Absorb the method
    spin_system.control.method=control.method;
    control=rmfield(control,'method');

    % Refuse to run Newton for trapezium integrator
    if strcmp(spin_system.control.method,'newton')&&...
       strcmp(spin_system.control.integrator,'trapezium')
        error('Newton optimiser is not available for trapezium integrator.');
    end

    % Refuse to run Goodwin for trapezium integrator
    if strcmp(spin_system.control.method,'goodwin')&&...
       strcmp(spin_system.control.integrator,'trapezium')
        error('Goodwin optimiser is not available for trapezium integrator.');
    end
    
else    

    % Default is LBFGS
    spin_system.control.method='lbfgs';
    
end

% Inform the user
report(spin_system,[pad('Optimisation method',60) spin_system.control.method]);

% Process iteration limit
if isfield(control,'max_iter')

    % Input validation
    if (~isnumeric(control.max_iter))||(~isreal(control.max_iter))||...
       (~isscalar(control.max_iter))||(mod(control.max_iter,1)~=0)||...
       (control.max_iter<0)
        error('control.max_iter must be a non-negative real integer.');
    end
    
    % Absorb the specification
    spin_system.control.max_iter=control.max_iter;
    control=rmfield(control,'max_iter');
    
else
    
    % Default is 100
    spin_system.control.max_iter=100;
    
end

% Inform the user
report(spin_system,[pad('Maximum number of iterations',60) ...
                    pad(num2str(spin_system.control.max_iter),20)]);

% Process step norm tolerance
if isfield(control,'tol_x')

    % Input validation
    if (~isnumeric(control.tol_x))||(~isreal(control.tol_x))||...
       (~isscalar(control.tol_x))||(control.tol_x<0)
        error('control.tol_x must be a non-negative real scalar.');
    end
    
    % Absorb the spec
    spin_system.control.tol_x=control.tol_x;
    control=rmfield(control,'tol_x');
    
else
    
    % Default is 0.1%
    spin_system.control.tol_x=1e-3;

end

% Inform the user
report(spin_system,[pad('Termination tolerance on |delta_x|',60) ...
                    pad(num2str(spin_system.control.tol_x,'%0.8g'),20)]);

% Process gradient norm tolerance
if isfield(control,'tol_g')

    % Input validation
    if (~isnumeric(control.tol_g))||(~isreal(control.tol_g))||...
       (~isscalar(control.tol_g))||(control.tol_g<0)
        error('control.tol_g must be a non-negative real scalar.');
    end
    
    % Absorb the spec
    spin_system.control.tol_g=control.tol_g;
    control=rmfield(control,'tol_g');
    
else
    
    % Default is 1e-6
    spin_system.control.tol_g=1e-6;
    
end

% Inform the user
report(spin_system,[pad('Termination tolerance on |grad(x)|',60) ...
                    pad(num2str(spin_system.control.tol_g,'%0.8g'),20)]);

% Set up LBFGS history
if strcmp(spin_system.control.method,'lbfgs')
    
    % Decide history length
    if isfield(control,'n_grads')

        % Input validation
        if (~isnumeric(control.n_grads))||(~isreal(control.n_grads))||...
           (~isscalar(control.n_grads))||(mod(control.n_grads,1)~=0)||...
           (control.n_grads<2)
            error('control.n_grads must be a real integer gerater than 1.');
        end
        
        % Absorb the spec
        spin_system.control.n_grads=control.n_grads;
        control=rmfield(control,'n_grads');
        
    else
        
        % Default is 10 gradients
        spin_system.control.n_grads=10;
        
    end
    
    % Inform the user
    report(spin_system,[pad('Number of gradients in LBFGS history',60) ...
                        pad(num2str(spin_system.control.n_grads),20)]);
    
end

% Process phase cycle
if isfield(control,'phase_cycle')

    % Input validation
    if (~isnumeric(control.phase_cycle))||(~isreal(control.phase_cycle))
        error('control.phase_cycle must be a real numeric array.');
    end
    if mod(numel(spin_system.control.operators),2)~=0
        error('phase cycles require an even number of control operators.');
    end
    if size(control.phase_cycle,2)~=numel(spin_system.control.operators)/2+2
        error(['each phase cycle line must have '...
               int2str(numel(spin_system.control.operators)/2+2) ' elements.']);
    end

    % Absorb phase cycle specification
    spin_system.control.phase_cycle=control.phase_cycle;
    control=rmfield(control,'phase_cycle');
    
    % Inform the user
    report(spin_system,[pad('Number of phase cycle elements',60) ...
                        int2str(size(spin_system.control.phase_cycle,1))]);
    
else
    
    % Empty array with one row
    spin_system.control.phase_cycle=zeros(1,0);
    
    % Inform the user
    report(spin_system,[pad('Number of phase cycle elements',60) '0']);
    
end

% Process basis set
if isfield(control,'basis')

    % Input validation
    if (~isnumeric(control.basis))||(~isreal(control.basis))
        error('control.basis must be a real matrix.');
    end
    if size(control.basis,2)~=spin_system.control.pulse_ntpts
        error(['the number of columns in control.basis must be ' ...
               int2str(spin_system.control.pulse_ntpts)]);
    end
    overlap_matrix=control.basis*control.basis';
    if norm(overlap_matrix-eye(size(overlap_matrix)),1)>1e-6
        error('waveform basis is not orthonormal.');
    end

    % Absorb waveform basis
    spin_system.control.basis=control.basis; 
    control=rmfield(control,'basis');

    % Inform the user
    report(spin_system,[pad('Number of waveform basis functions',60) ...
                        int2str(size(spin_system.control.basis,1))]);
                    
else

    % Empty array for the basis set
    spin_system.control.basis=[];

end

% Process keyhole schedule
if isfield(control,'keyholes')

    % Input validation
    if strcmp(spin_system.control.integrator,'trapezium')
        error('keyholes are not supported with trapezium integrator.');
    end
    if ~iscell(control.keyholes)
        error('control.keyholes must be a cell array.');
    end
    for n=1:numel(control.keyholes)
        if (~isempty(control.keyholes{n}))&&(~isa(control.keyholes{n},'function_handle'))
            error('non-empty elements of control.keyholes must be function handles.');
        end
        if isa(control.keyholes{n},'function_handle')
            a=randn(1,size(spin_system.control.rho_init{1},1));
            b=randn(1,size(spin_system.control.rho_init{1},1));
            keyhole=control.keyholes{n};
            if norm(keyhole(a)-keyhole(keyhole(a)),2)/norm(a,2)>1e-6
                error('keyhole functions must be idempotent.');
            end
            if norm(keyhole(a)+keyhole(b)-keyhole(a+b),2)/norm(a+b,2)>1e-6
                error('keyhole functions must be linear.');
            end
        end
    end
    
    % Absorb keyhole schedule
    spin_system.control.keyholes=control.keyholes;
    control=rmfield(control,'keyholes');
    
    % Inform the user
    report(spin_system,[pad('Keyholes present',60) 'yes']);
    
else
    
    % Set up empty keyhole schedule
    spin_system.control.keyholes=cell(1,spin_system.control.pulse_ntpts);
    
    % Inform the user
    report(spin_system,[pad('Keyholes present',60) 'no']);
    
end

% Process penalties                
if isfield(control,'penalties')&&isfield(control,'p_weights')

    % Input validation
    if (~iscell(control.penalties))||any(~cellfun(@ischar,control.penalties(:)))
        error('control.penalties must be a cell array of character strings.');
    end
    if any(~ismember(control.penalties(:),{'none','NS','SNS','DNS','SNSA'}))
        error('elements of control.penalties can be ''none'',''NS'', ''SNS'', ''SNSA'', or ''DNS''');
    end
    if (~isnumeric(control.p_weights))||(~isreal(control.p_weights))||...
       (~isrow(control.p_weights))||any(control.p_weights(:)<0)
            error('control.p_weights must be a row vector of non-negative real numbers.');
    end
    
    % Absorb penalties and their weights
    spin_system.control.penalties=control.penalties; control=rmfield(control,'penalties');
    spin_system.control.p_weights=control.p_weights; control=rmfield(control,'p_weights');
    
else
    
    % Default is SNS penalty
    spin_system.control.penalties={'SNS'};
    spin_system.control.p_weights=100.0;
    
end

% Inform the user
for n=1:numel(spin_system.control.penalties)
    report(spin_system,[pad(['Penalty function ' num2str(n) ', weight ' ...
                              num2str(spin_system.control.p_weights(n),'%.9g')],60)...
                              spin_system.control.penalties{n}]);
end

% Process lower and upper bounds
if isfield(control,'u_bound')&&isfield(control,'l_bound')

    % Input validation
    if (~isnumeric(control.u_bound))||(~isreal(control.u_bound))||(~isscalar(control.u_bound))
        error('control.u_bound must be a real scalar.');
    end
    if (~isnumeric(control.l_bound))||(~isreal(control.l_bound))||(~isscalar(control.l_bound))
        error('control.l_bound must be a real scalar.');
    end
        
    % Absorb bounds
    spin_system.control.u_bound=control.u_bound; control=rmfield(control,'u_bound');
    spin_system.control.l_bound=control.l_bound; control=rmfield(control,'l_bound');  
    
else
    
    % Default is nominal power
    spin_system.control.u_bound=+1;
    spin_system.control.l_bound=-1;  
    
end

% Inform the user
report(spin_system,[pad('SNS penalty ceiling, fraction of power level',60) ...
                         num2str(spin_system.control.u_bound,'%+.9g')]);
report(spin_system,[pad('SNS penalty floor, fraction of power level',60) ...
                         num2str(spin_system.control.l_bound,'%+.9g')]);

% Process ensemble correlations
if isfield(control,'ens_corrs')

    % Input validation (more inside ensemble function)
    if (~iscell(control.ens_corrs))||(~all(cellfun(@ischar,control.ens_corrs)))
        error('control.ens_corrs must be a cell array of character strings.');
    end
    for n=1:numel(control.ens_corrs)
        if ~ismember(control.ens_corrs{n},{'rho_ens','rho_drift','power_drift'})
            error('control.ens_corrs contains an unknown ensemble correlation setting')
        end
    end
    
    % Absorb ensemble correlations
    spin_system.control.ens_corrs=control.ens_corrs;
    control=rmfield(control,'ens_corrs');

else

    % No ensemble correlations
    spin_system.control.ens_corrs={};

end

% Inform the user
for n=1:numel(spin_system.control.ens_corrs)
    report(spin_system,[pad('Ensemble correlation',60) ...
                        spin_system.control.ens_corrs{n}]);
end

% Parallelisation strategy
if isfield(control,'parallel')

    % Input validation
    if ~ischar(control.parallel)
        error('control.parallel must be a character string.');
    end
    if ~ismember(control.parallel,{'ensemble','time'})
        error('control.parallel must be ''ensemble'' or ''time''.');
    end

    % Absorb the specification
    spin_system.control.parallel=control.parallel;
    control=rmfield(control,'parallel');

else

    % Default is over time steps
    spin_system.control.parallel='time';

end

% Inform the user
report(spin_system,[pad('Parallelisation strategy',60) spin_system.control.parallel]);
                     
% Plotting options
if isfield(control,'plotting')

    % Input validation
    if (~iscell(control.plotting))||any(~cellfun(@ischar,control.plotting(:)))
        error('control.plotting must be a cell array of character strings.');
    end
    if any(~ismember(control.plotting(:),{'xy_controls','phi_controls',...
                                          'amp_controls','correlation_order',...
                                          'coherence_order','local_each_spin',...
                                          'total_each_spin','level_populations',...
                                          'robustness','spectrogram','time_by_slice'}))
        error('unrecognised element encountered in control.plotting');
    end

    % Absorb the specification
    spin_system.control.plotting=control.plotting;
    control=rmfield(control,'plotting');
    
else
    
    % Default is no plotting
    spin_system.control.plotting={};
    
end

% Distortion before plotting
if isfield(control,'distplot')

    % Input validation
    if (~iscell(control.distplot))||...
        any(~cellfun(@(x)isa(x,'function_handle'),control.distplot(:)))
        error('control.distplot must be a cell array of function handles.');
    end
    
    % Absorb the specification
    spin_system.control.distplot=control.distplot;
    control=rmfield(control,'distplot');
    
else
    
    % Plot ideal waveform by default
    spin_system.control.distplot={};
    
end

% Process trajectory options
if isfield(control,'traj_opts')

    % Input validation
    if (~iscell(control.traj_opts))||any(~cellfun(@ischar,control.traj_opts(:)))
        error('control.traj_opts must be a cell array of character strings.');
    end
    if any(~ismember(control.traj_opts(:),{'average'}))
        error('unrecognised element encountered in control.traj_opts');
    end
   
    % Absorb the specification
    spin_system.control.traj_opts=control.traj_opts;
    control=rmfield(control,'traj_opts');
    
else
    
    % Default is no sum
    spin_system.control.traj_opts={};
    
end

% Process checkpoint file
if isfield(control,'checkpoint')

    % Input validation
    if ~ischar(control.checkpoint)
        error('control.checkpoint must be a character string.');
    end
    
    % Store checkpoint file name
    spin_system.control.checkpoint=control.checkpoint;
    control=rmfield(control,'checkpoint');
    
    % Inform the user
    report(spin_system,[pad('Checkpoint file name (scratch directory)',60) ...
                        spin_system.control.checkpoint]);
                    
end

% Process video file name
if (~isempty(spin_system.control.plotting))&&...
     isfield(control,'video_file')

    % Input validation
    if ~ischar(control.video_file)
        error('control.video_file must be a character string.');
    end

    % Absorb the video file name
    spin_system.control.video_file=control.video_file;
    control=rmfield(control,'video_file');
    
    % Inform the user
    report(spin_system,[pad('Video file name',60) ...
                        spin_system.control.video_file]);
                 
end

% In Liouville space, warn about non-Hermitian controls
if ~all(cellfun(@ishermitian,spin_system.control.operators))
    report(spin_system,'WARNING: not all control operators are Hermitian');
end

% Warn that plotting is expensive
if ~isempty(spin_system.control.plotting)
    report(spin_system,'WARNING: diagnostic plots are computationally expensive');
end

% Warn that StStSt ignores the initial condition
if spin_system.control.steady
    report(spin_system,'WARNING: stroboscopic steady state mode, rho_init ignored.');
end

% Internal parameters that users do not need to know about
spin_system.control.ls_c1=1e-2;        % Line search: sufficient decrease condition
spin_system.control.ls_c2=0.9;         % Line search: curvature condition on grad
spin_system.control.ls_tau1=3;         % Line search: bracket expansion factor
spin_system.control.ls_tau2=0.1;       % Line search: left section contraction
spin_system.control.ls_tau3=0.5;       % Line search: right section contraction
spin_system.control.reg_max_iter=2500; % RFO: max regularisation iterations
spin_system.control.reg_alpha=1;       % RFO: scaling factor
spin_system.control.reg_phi=0.9;       % RFO: conditioning multiplier
spin_system.control.reg_max_cond=1e4;  % RFO: max condition number

% Accept pulse sequence parameters
if isfield(control,'parameters')
    spin_system.control.parameters=control.parameters;
    control=rmfield(control,'parameters');
    report(spin_system,'control.parameters received without parsing.');
end

% Catch unparsed fields
unparsed=fieldnames(control);
if ~isempty(unparsed)
    for n=1:numel(unparsed)
        report(spin_system,['ERROR: unrecognised or mismatched option - ' unparsed{n}]);
    end
    error('there are problems with the control structure.');
end

end

% Consistency enforcement
function grumble(spin_system,control) %#ok<INUSD> 

% Add further sanity checks here

end

% Whenever anyone accuses some person of being 'unfeeling' he means that
% that person is just. He means that that person has no causeless emotions
% and will not grant him a feeling which he does not deserve [...] justice
% is the opposite of charity.
%
% Ayn Rand, "Atlas Shrugged"

