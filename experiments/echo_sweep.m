% Two-pulse echo-detected frequency-swept experiment, static or under
% magic angle spinning, in Hilbert space, written for the EPR case of
% a spinning P1 centre in diamond. Two pulses of equal duration are
% separated by a delay, the carrier is stepped across the sweep, and
% the complex echo integral is returned at each carrier offset. The sequence steps through the Hamiltonian rotor stack
% supplied by singlerot.m: at each time step, the stack element near-
% est to the rotor phase at the middle of the step is used, and the
% rotor phase at the start of the sequence, which stands in for the
% crystallite azimuth about the rotor axis, is averaged over. The co-
% herence pathway of the pulsed spin (-1 after the first pulse, +1
% after the second) is selected in place of a phase cycle. Syntax:
%
%           echo=echo_sweep(spin_system,parameters,H,R,K)
%
% Parameters:
%
%    parameters.spins     - one-element cell array naming the spin
%                           the pulses are applied to, e.g. {'E'}
%
%    parameters.rho0      - initial state, a density matrix
%
%    parameters.coil      - detection state, a density matrix
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
%    parameters.timestep  - propagation time step, seconds; the
%                           pulses, the delay, and the echo win-
%                           dow are rounded to whole steps
%
%    parameters.rate      - spinning rate, Hz, zero for a static
%                           sample
%
%    parameters.nphases   - number of rotor phases at the start
%                           of the sequence to average over
%
%    parameters.sweep     - width of the carrier sweep, Hz
%
%    parameters.npoints   - number of carrier offsets, placed on
%                           the ft_axis grid of the sweep
%
%    parameters.spc_dim   - number of elements in the rotor stack,
%                           received from context function
%
%    H  - vector cell array of Hamiltonian matrices, one for each
%         rotor phase, received from context function
%
%    R  - relaxation superoperator, received from context func-
%         tion, not used
%
%    K  - kinetics superoperator, received from context function,
%         not used
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

function echo=echo_sweep(spin_system,parameters,H,~,~)

% Check consistency
grumble(spin_system,parameters,H);

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

% Rotor stack advance at the middle of each time step
stack_shift=round(parameters.rate*parameters.timestep*((1:nsteps)-1/2)*parameters.spc_dim);

% Rotor stack indices at the start of the sequence
start_idx=floor((0:(parameters.nphases-1))*parameters.spc_dim/parameters.nphases);

% Free evolution propagators at every rotor phase
p_free=cell(parameters.spc_dim,1);
for n=1:parameters.spc_dim
    p_free{n}=propagator(spin_system,H{n},parameters.timestep);
end

% Preallocate the answer
echo=zeros(parameters.npoints,1);

% Loop over carrier offsets
for k=1:parameters.npoints

    % Carrier offset propagator
    p_off=propagator(spin_system,2*pi*offsets(k)*sz,parameters.timestep);

    % Pulse propagators at every rotor phase
    p_pulse=cell(parameters.spc_dim,1);
    for n=1:parameters.spc_dim
        p_pulse{n}=propagator(spin_system,H{n}+2*pi*offsets(k)*sz+...
                              2*pi*parameters.pulse_frq*sx,parameters.timestep);
    end

    % Loop over rotor phases at the start of the sequence
    for j=1:parameters.nphases

        % Rotor stack indices at each time step
        idx=mod(start_idx(j)+stack_shift,parameters.spc_dim)+1;

        % First pulse
        rho=parameters.rho0;
        for s=1:pulse_steps
            rho=p_pulse{idx(s)}*rho*p_pulse{idx(s)}';
        end

        % Select the -1 coherence on the pulsed spin
        rho=coherence(spin_system,rho,{{parameters.spins{1},-1}});

        % Interpulse delay
        for s=(pulse_steps+1):(pulse_steps+delay_steps)
            rho=p_off*p_free{idx(s)}*rho*p_free{idx(s)}'*p_off';
        end

        % Second pulse
        for s=(pulse_steps+delay_steps+1):(2*pulse_steps+delay_steps)
            rho=p_pulse{idx(s)}*rho*p_pulse{idx(s)}';
        end

        % Select the +1 coherence on the pulsed spin
        rho=coherence(spin_system,rho,{{parameters.spins{1},+1}});

        % Integrate the signal over the echo window
        for s=(2*pulse_steps+delay_steps+1):nsteps
            rho=p_off*p_free{idx(s)}*rho*p_free{idx(s)}'*p_off';
            echo(k)=echo(k)+trace(parameters.coil'*rho);
        end

    end

end

% Integrate over the echo window and average over the rotor phases
echo=parameters.timestep*echo/parameters.nphases;

end

% Consistency enforcement
function grumble(spin_system,parameters,H)
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


