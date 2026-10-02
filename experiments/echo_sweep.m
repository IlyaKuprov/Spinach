% Two-pulse echo-detected frequency-swept experiment, static or under
% magic angle spinning, in Hilbert or Liouville space, for the EPR case of
% a spinning P1 centre in diamond. Two pulses of equal duration are
% separated by a delay, the carrier is stepped across the sweep, and
% the complex echo integral is returned at each carrier offset.
% In Hilbert space the sequence steps through the Hamiltonian stack
% supplied by singlerot.m: at each time step, the stack element nearest
% to the rotor phase at the middle of the step is used, the phase de-
% creasing with time for a positive rate as in the Liouville space
% branch of singlerot.m, and only the elements visited are exponen-
% tiated. The rotor phase at the start of the sequence, which stands
% in for the crystallite azimuth about the rotor axis, is averaged
% over. In Liouville space, singlerot.m builds the Fokker-Planck
% generator and handles powder rotor-phase averaging; only the echo
% detection window is sampled. The coherence pathway (-1 after the first
% pulse, +1 after the second) is selected in place of a phase cycle.
% Syntax:
%
%           echo=echo_sweep(spin_system,parameters,H,R,K)
%
% Parameters:
%
%    parameters.spins     - one-element cell array naming the spin
%                           the pulses are applied to, e.g. {'E'}
%
%    parameters.rho0      - initial density matrix (Hilbert) or state
%                           vector (Liouville)
%
%    parameters.coil      - detection matrix (Hilbert) or vector
%                           (Liouville)
%
%    parameters.pulse_dur - duration of each pulse, seconds
%
%    parameters.pulse_frq - nutation frequency of the pulses, Hz
%
%    parameters.tau       - delay between the end of the first
%                           pulse and the start of the second
%                           pulse, seconds
%
%    parameters.echo_win  - echo integration window after the end
%                           of the second pulse, seconds
%
%    parameters.timestep  - time step, seconds; in Hilbert space
%                           pulses, delay, and window are rounded
%                           to steps; in Liouville space only the
%                           echo window is sampled
%
%    parameters.rate      - spinning rate, Hz, zero for a static
%                           sample
%
%    parameters.nphases   - Hilbert only: number of initial rotor
%                           phases to average over
%
%    parameters.sweep     - width of the carrier sweep, Hz
%
%    parameters.npoints   - number of carrier offsets, placed on
%                           the ft_axis grid of the sweep
%
%    parameters.spc_dim   - number of rotor grid points supplied
%                           by singlerot.m
%
%    H  - Hilbert: cell array of rotor-phase Hamiltonians; Liou-
%         ville: rotor-augmented Fokker-Planck generator
%
%    R  - relaxation superoperator, used in Liouville space
%
%    K  - kinetics superoperator, used in Liouville space
%
% Outputs:
%
%    echo - complex echo signal integrated over the echo window (a
%           sum over the time steps multiplied by the time step) and
%           averaged over the rotor phases at the start of the sequ-
%           ence, at each carrier offset, a column vector with
%           parameters.npoints elements
%
% Note: the elements of the rotor stack must commute with the Lz op-
%       erator of the pulsed spin, as they do for an electron under
%       the 'esr' assumption set, because the carrier offset is app-
%       lied as a separate propagator and only the pulse propagators
%       are rebuilt at each carrier offset.
%
% Note: the rotor stack is a table of the Hamiltonian against the
%       rotor phase, its resolution should match the time step at
%       the fastest spinning rate used, parameters.max_rank of the
%       context function of about 1/(2*abs(rate)*timestep) there;
%       a finer stack costs propagators without gaining accuracy be-
%       yond the time step, a coarser one loses rotor phase resolu-
%       tion. At slower rates, consecutive steps reuse elements.
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=echo_sweep.m>

function echo=echo_sweep(spin_system,parameters,H,R,K)

% Check consistency
grumble(spin_system,parameters,H,R,K);

% Fokker-Planck singlerot supplies a rotor-augmented Liouvillian
if ismember(spin_system.bas.formalism,{'sphten-liouv','zeeman-liouv'})

    % Carrier offsets and microwave operators in the rotor-augmented space
    offsets=ft_axis(0,parameters.sweep,parameters.npoints);
    sx=kron(speye(parameters.spc_dim),operator(spin_system,'Lx',parameters.spins{1}));
    sz=kron(speye(parameters.spc_dim),operator(spin_system,'Lz',parameters.spins{1}));
    L=H+1i*R+1i*K;

    % Keep the detection state and accumulator on the propagator device
    coil=parameters.coil;
    echo=zeros(parameters.npoints,1);
    if ismember('gpu',spin_system.sys.enable)
        coil=gpuArray(coil);
        echo=gpuArray(echo);
    end
    echo_steps=round(parameters.echo_win/parameters.timestep);
    for k=1:parameters.npoints

        % The free generator and the finite-pulse generator at this carrier
        L0=L+2*pi*offsets(k)*sz;
        Lp=L0+2*pi*parameters.pulse_frq*sx;

        % First pulse, coherence selection, and interpulse delay
        rho=step(spin_system,Lp,parameters.rho0,parameters.pulse_dur);
        rho=coherence(spin_system,rho,{{parameters.spins{1},-1}});
        if parameters.tau>0

            % Match the CPU or GPU norm used by step for long delays
            if ismember('gpu',spin_system.sys.enable)&&~isa(L0,'polyadic')
                generator_norm=norm(L0,inf);
            else
                generator_norm=cheap_norm(L0);
            end
            nsteps=max(1,ceil(generator_norm*parameters.tau/2e4));
            for n=1:nsteps
                rho=step(spin_system,L0,rho,parameters.tau/nsteps);
            end
        end

        % Second pulse and refocused coherence
        rho=step(spin_system,Lp,rho,parameters.pulse_dur);
        rho=coherence(spin_system,rho,{{parameters.spins{1},+1}});

        % Only the finite detection window requires time samples
        for n=1:echo_steps
            rho=step(spin_system,L0,rho,parameters.timestep);
            echo(k)=echo(k)+coil'*rho;
        end
    end

    % Echo integral has units of signal times seconds
    echo=echo*parameters.timestep;

    % Return the completed spectrum to CPU memory for powder averaging
    if isa(echo,'gpuArray')
        echo=gather(echo);
    end

    return
end

% Carrier offsets across the sweep
offsets=ft_axis(0,parameters.sweep,parameters.npoints);

% Pulse and offset operators
sx=operator(spin_system,'Lx',parameters.spins{1});
sz=operator(spin_system,'Lz',parameters.spins{1});

% Step counts of the pulses, the delay, and the echo window
pulse_steps=round(parameters.pulse_dur/parameters.timestep);
delay_steps=round(parameters.tau/parameters.timestep);
echo_steps=round(parameters.echo_win/parameters.timestep);
nsteps=2*pulse_steps+delay_steps+echo_steps;

% Rotor stack shift at the middle of each time step
stack_shift=round(parameters.rate*parameters.timestep*((1:nsteps)-1/2)*parameters.spc_dim);

% Rotor stack indices at the start of the sequence
start_idx=floor((0:(parameters.nphases-1))*parameters.spc_dim/parameters.nphases);

% Rotor stack indices at each time step for each start phase
idx=mod(start_idx'-stack_shift,parameters.spc_dim)+1;

% Rotor stack elements visited by the sequence
used=unique(idx(:))';

% Free evolution propagators at the rotor phases visited
p_free=cell(parameters.spc_dim,1);
for n=used
    p_free{n}=propagator(spin_system,H{n},parameters.timestep);
end

% Preallocate the answer
echo=zeros(parameters.npoints,1);

% Loop over carrier offsets
for k=1:parameters.npoints

    % Carrier offset propagator
    p_off=propagator(spin_system,2*pi*offsets(k)*sz,parameters.timestep);

    % Pulse propagators at the rotor phases visited
    p_pulse=cell(parameters.spc_dim,1);
    for n=used
        p_pulse{n}=propagator(spin_system,H{n}+2*pi*offsets(k)*sz+...
                              2*pi*parameters.pulse_frq*sx,parameters.timestep);
    end

    % Loop over rotor phases at the start of the sequence
    for j=1:parameters.nphases

        % First pulse
        rho=parameters.rho0;
        for s=1:pulse_steps
            rho=p_pulse{idx(j,s)}*rho*p_pulse{idx(j,s)}';
        end

        % Select the -1 coherence on the pulsed spin
        rho=coherence(spin_system,rho,{{parameters.spins{1},-1}});

        % Interpulse delay
        for s=(pulse_steps+1):(pulse_steps+delay_steps)
            rho=p_off*p_free{idx(j,s)}*rho*p_free{idx(j,s)}'*p_off';
        end

        % Second pulse
        for s=(pulse_steps+delay_steps+1):(2*pulse_steps+delay_steps)
            rho=p_pulse{idx(j,s)}*rho*p_pulse{idx(j,s)}';
        end

        % Select the +1 coherence on the pulsed spin
        rho=coherence(spin_system,rho,{{parameters.spins{1},+1}});

        % Integrate the signal over the echo window
        for s=(2*pulse_steps+delay_steps+1):nsteps
            rho=p_off*p_free{idx(j,s)}*rho*p_free{idx(j,s)}'*p_off';
            echo(k)=echo(k)+trace(parameters.coil'*rho);
        end

    end

end

% Integrate over the echo window and average over the rotor phases
echo=parameters.timestep*echo/parameters.nphases;

end

% Consistency enforcement
function grumble(spin_system,parameters,H,R,K)

if ismember(spin_system.bas.formalism,{'sphten-liouv','zeeman-liouv'})
    if (~isnumeric(H))||(~isnumeric(R))||(~isnumeric(K))||...
       (~ismatrix(H))||(~isequal(size(H),size(R),size(K)))||...
       (size(H,1)~=size(H,2))||...
       (~allfinite(H))||(~allfinite(R))||...
       (~allfinite(K))
        error('H, R, and K must be finite, square, equal-sized matrices.');
    end

    required={'spc_dim','spins','rho0','coil','pulse_dur','pulse_frq',...
              'tau','echo_win','timestep','sweep','npoints'};
    for n=1:numel(required)
        if ~isfield(parameters,required{n})
            error('parameters.%s is required for echo_sweep.',required{n});
        end
    end
    validateattributes(parameters.spc_dim,{'numeric'},...
                       {'scalar','real','finite','integer','positive'});
    if (~iscell(parameters.spins))||(numel(parameters.spins)~=1)||...
       (~ischar(parameters.spins{1}))||isempty(parameters.spins{1})
        error('parameters.spins must contain one pulsed isotope.');
    end
    if mod(size(H,1),parameters.spc_dim)~=0
        error('Liouvillian dimension must be divisible by parameters.spc_dim.');
    end
    for name={'rho0','coil'}
        state_vec=parameters.(name{1});
        if (~isnumeric(state_vec))||(~iscolumn(state_vec))||...
           (numel(state_vec)~=size(H,1))||(~all(isfinite(nonzeros(state_vec))))
            error('parameters.%s must be a finite rotor-augmented state.',name{1});
        end
    end
    for name={'pulse_dur','pulse_frq','echo_win','timestep','sweep'}
        validateattributes(parameters.(name{1}),{'numeric'},...
                           {'scalar','real','finite','positive'});
    end
    validateattributes(parameters.tau,{'numeric'},...
                       {'scalar','real','finite','nonnegative'});
    validateattributes(parameters.npoints,{'numeric'},...
                       {'scalar','real','finite','integer','>',2});
    if round(parameters.echo_win/parameters.timestep)<1
        error('parameters.echo_win must contain at least one time sample.');
    end

    return
end
if ~strcmp(spin_system.bas.formalism,'zeeman-hilb')
    error('this function is only available in zeeman-hilb formalism.');
end
if ~isfield(parameters,'spc_dim')
    error('rotor stack size must be specified in parameters.spc_dim field.');
end
if (~isnumeric(parameters.spc_dim))||(~isreal(parameters.spc_dim))||...
   (~isscalar(parameters.spc_dim))||(mod(parameters.spc_dim,1)~=0)||...
   (parameters.spc_dim<1)
    error('parameters.spc_dim must be a positive real integer.');
end
if (~iscell(H))||(~isvector(H))||(numel(H)~=parameters.spc_dim)||...
   (~all(cellfun(@(x)isnumeric(x)&&ismatrix(x)&&(size(x,1)==size(x,2)),H)))
    error('H must be a vector cell array of parameters.spc_dim square matrices.');
end
if ~all(cellfun(@(x)all(size(x)==size(H{1})),H))
    error('all matrices in H must have the same dimension.');
end
if ~all(cellfun(@(x)all(isfinite(nonzeros(x))),H))
    error('the elements of H must have finite entries.');
end
if any(cellfun(@(x)norm(x-x',1)>spin_system.tols.liouv_zero*norm(x,1),H))
    error('the elements of H must be Hermitian.');
end
if ~isfield(parameters,'spins')
    error('the pulsed spin must be specified in parameters.spins field.');
end
if (~iscell(parameters.spins))||(numel(parameters.spins)~=1)||...
   (~ischar(parameters.spins{1}))
    error('parameters.spins must be a one-element cell array of character strings.');
end
sz=operator(spin_system,'Lz',parameters.spins{1});
if ~isequal(size(H{1}),size(sz))
    error('the elements of H must have the dimension of the spin system.');
end
if any(cellfun(@(x)norm(x*sz-sz*x,1)>spin_system.tols.liouv_zero,H))
    error('the elements of H must commute with the Lz operator of the pulsed spin.');
end
if ~isfield(parameters,'rho0')
    error('initial state must be specified in parameters.rho0 field.');
end
if (~isnumeric(parameters.rho0))||(~isequal(size(parameters.rho0),size(H{1})))||...
   (~all(isfinite(nonzeros(parameters.rho0))))
    error('parameters.rho0 must be a finite matrix of the same dimension as the elements of H.');
end
if ~isfield(parameters,'coil')
    error('detection state must be specified in parameters.coil field.');
end
if (~isnumeric(parameters.coil))||(~isequal(size(parameters.coil),size(H{1})))||...
   (~all(isfinite(nonzeros(parameters.coil))))
    error('parameters.coil must be a finite matrix of the same dimension as the elements of H.');
end
if ~isfield(parameters,'pulse_dur')
    error('pulse duration must be specified in parameters.pulse_dur field.');
end
if (~isnumeric(parameters.pulse_dur))||(~isreal(parameters.pulse_dur))||...
   (~isscalar(parameters.pulse_dur))||(~isfinite(parameters.pulse_dur))||(parameters.pulse_dur<=0)
    error('parameters.pulse_dur must be a finite positive real scalar.');
end
if ~isfield(parameters,'pulse_frq')
    error('pulse nutation frequency must be specified in parameters.pulse_frq field.');
end
if (~isnumeric(parameters.pulse_frq))||(~isreal(parameters.pulse_frq))||...
   (~isscalar(parameters.pulse_frq))||(~isfinite(parameters.pulse_frq))||(parameters.pulse_frq<=0)
    error('parameters.pulse_frq must be a finite positive real scalar.');
end
if ~isfield(parameters,'tau')
    error('interpulse delay must be specified in parameters.tau field.');
end
if (~isnumeric(parameters.tau))||(~isreal(parameters.tau))||...
   (~isscalar(parameters.tau))||(~isfinite(parameters.tau))||(parameters.tau<0)
    error('parameters.tau must be a finite non-negative real scalar.');
end
if ~isfield(parameters,'echo_win')
    error('echo integration window must be specified in parameters.echo_win field.');
end
if (~isnumeric(parameters.echo_win))||(~isreal(parameters.echo_win))||...
   (~isscalar(parameters.echo_win))||(~isfinite(parameters.echo_win))||(parameters.echo_win<=0)
    error('parameters.echo_win must be a finite positive real scalar.');
end
if ~isfield(parameters,'timestep')
    error('propagation time step must be specified in parameters.timestep field.');
end
if (~isnumeric(parameters.timestep))||(~isreal(parameters.timestep))||...
   (~isscalar(parameters.timestep))||(~isfinite(parameters.timestep))||(parameters.timestep<=0)
    error('parameters.timestep must be a finite positive real scalar.');
end
if parameters.pulse_dur<parameters.timestep
    error('parameters.pulse_dur must not be shorter than parameters.timestep.');
end
if parameters.echo_win<parameters.timestep
    error('parameters.echo_win must not be shorter than parameters.timestep.');
end
if ~isfield(parameters,'rate')
    error('spinning rate must be specified in parameters.rate field.');
end
if (~isnumeric(parameters.rate))||(~isreal(parameters.rate))||...
   (~isscalar(parameters.rate))||(~isfinite(parameters.rate))
    error('parameters.rate must be a finite real scalar.');
end
if ~isfield(parameters,'nphases')
    error('number of rotor phases must be specified in parameters.nphases field.');
end
if (~isnumeric(parameters.nphases))||(~isreal(parameters.nphases))||...
   (~isscalar(parameters.nphases))||(mod(parameters.nphases,1)~=0)||...
   (parameters.nphases<1)||(parameters.nphases>parameters.spc_dim)
    error('parameters.nphases must be a positive real integer not exceeding parameters.spc_dim.');
end
if ~isfield(parameters,'sweep')
    error('carrier sweep width must be specified in parameters.sweep field.');
end
if (~isnumeric(parameters.sweep))||(~isreal(parameters.sweep))||...
   (~isscalar(parameters.sweep))||(~isfinite(parameters.sweep))||(parameters.sweep<=0)
    error('parameters.sweep must be a finite positive real scalar.');
end
if ~isfield(parameters,'npoints')
    error('number of carrier offsets must be specified in parameters.npoints field.');
end
if (~isnumeric(parameters.npoints))||(~isreal(parameters.npoints))||...
   (~isscalar(parameters.npoints))||(mod(parameters.npoints,1)~=0)||...
   (parameters.npoints<=2)
    error('parameters.npoints must be a real integer greater than 2.');
end
end


