% Overtone DANTE experiment with frequency-domain acquisition.
%
%     spectrum=overtone_dante(spin_system,parameters,H,R,K)
%
% where H is the Hamiltonian matrix, R is the relaxation matrix
% and K is the chemical kinetics matrix. The following parameters
% are required:
%
%     parameters.pulse_dur - duration of the pulse, seconds
%
%     parameters.pulse_amp - amplitude of the pulse, rad/s
%
%     parameters.pulse_num - number of pulses within the 
%                            rotor period
%
%     parameters.n_periods - number of rotor periods that the 
%                            sequence is active for
%
%     parameters.spins     - overtone-active nucleus, specified
%                            as a single-element cell array
%
%     parameters.spc_dim   - Fokker-Planck spatial dimension
%
%     parameters.Lx        - X Zeeman operator on the
%                            quadrupolar nucleus
%
%     parameters.rf_frq    - pulse frequency offset from the
%                            overtone frequency, Hz
%
%     parameters.rate      - rotor frequency in Hz
%
%     parameters.sweep     - acquisition sweep range, Hz
%
%     parameters.npoints   - number of acquisition points
%
%     parameters.rho0      - initial condition, usually Lz
%
%     parameters.coil      - detection state, usually L+
%
% Outputs:
%
%     spectrum            - overtone spectrum
%
% ilya.kuprov@weizmann.ac.il
% m.carravetta@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=Overtone_dante.m>

function spectrum=overtone_dante(spin_system,parameters,H,R,K)

% Check consistency
grumble(parameters,H,R,K);

% Get the overtone frequency
ovt_frq=-2*spin(parameters.spins{1})*spin_system.inter.magnet/(2*pi);

% Project pulse operators
Lx=kron(speye(parameters.spc_dim),parameters.Lx);

% Timing parameters
rotor_period=abs(1/parameters.rate);
cycle_length=rotor_period/parameters.pulse_num;

% Bomb out if the schedule makes no sense
if (cycle_length-parameters.pulse_dur)<0
    error('DANTE pulse schedule does not fit into the rotor period.');
end

% Get the pulse frequency
omega=2*pi*ovt_frq-2*pi*parameters.rf_frq;

% Get the pulse Hamiltonian
pulseop=parameters.pulse_amp*Lx;
pulseop=average(spin_system,pulseop/2,H,pulseop/2,omega,'matrix_log');

% Precompute pulse propagator
P_pulse=propagator(spin_system,pulseop,parameters.pulse_dur);

% Precompute evolution propagator
P_evol=propagator(spin_system,H,cycle_length-parameters.pulse_dur);

% Loop over rotor cycles
for k=1:parameters.n_periods
    
    % Loop over pulses within the rotor cycle
    for n=1:parameters.pulse_num
        
        % Pulse and evolution
        parameters.rho0=P_evol*P_pulse*parameters.rho0;
        
    end
    
end
    
% Call the acquisition
spectrum=overtone_a(spin_system,parameters,H,R,K);

end

% Consistency enforcement
function grumble(parameters,H,R,K)
if (~isnumeric(H))||(~isnumeric(R))||(~isnumeric(K))||...
   (~ismatrix(H))||(~ismatrix(R))||(~ismatrix(K))
    error('H, R and K arguments must be matrices.');
end
if (~all(size(H)==size(R)))||(~all(size(R)==size(K)))
    error('H, R and K matrices must have the same dimension.');
end
if ~isfield(parameters,'spins')
    error('overtone-active nucleus must be specified in parameters.spins variable.');
elseif (~iscell(parameters.spins))||(numel(parameters.spins)~=1)
    error('parameters.spins should be a single-element cell array.');
end
if ~isfield(parameters,'spc_dim')
    error('Fokker-Planck dimension should be specified in parameters.spc_dim variable.');
elseif (~isnumeric(parameters.spc_dim))||(~isreal(parameters.spc_dim))||...
       (numel(parameters.spc_dim)~=1)||(parameters.spc_dim<1)||...
       (mod(parameters.spc_dim,1)~=0)
    error('parameters.spc_dim should be a positive integer.');
end
if ~isfield(parameters,'Lx')
    error('quadrupolar nucleus X pulse operator must be supplied in parameters.Lx field.');
end
if ~isfield(parameters,'rf_frq')
    error('parameters.rf_frq field is missing.');
elseif (~isnumeric(parameters.rf_frq))||(~isreal(parameters.rf_frq))||...
       (numel(parameters.rf_frq)~=1)
    error('parameters.rf_frq should be a real scalar.');
end
if ~isfield(parameters,'rate')
    error('rotor frequency must be specified in parameters.rate variable.');
elseif (~isnumeric(parameters.rate))||(~isreal(parameters.rate))||...
       (numel(parameters.rate)~=1)||(parameters.rate==0)
    error('parameters.rate should be a non-zero real scalar.');
end
if ~isfield(parameters,'npoints')
    error('number of points should be specified in parameters.npoints variable.');
end
if (~isnumeric(parameters.npoints))||(numel(parameters.npoints)~=1)||...
   (~isreal(parameters.npoints))||(parameters.npoints<1)||...
   (mod(parameters.npoints,1)~=0)
    error('parameters.npoints should be a positive integer.');
end
if ~isfield(parameters,'rho0')
    error('initial state must be specified in parameters.rho0 variable.');
end
if ~isfield(parameters,'coil')
    error('detection state must be specified in parameters.coil variable.');
end
if ~isfield(parameters,'sweep')
    error('spectral range must be specified in parameters.sweep variable.');
end
if (~isnumeric(parameters.sweep))||(~isreal(parameters.sweep))||...
   (numel(parameters.sweep)~=2)
    error('parameters.sweep vector must have two real elements.');
end
if ~isfield(parameters,'pulse_dur')
    error('pulse duration must be supplied in parameters.pulse_dur field.');
elseif (~isnumeric(parameters.pulse_dur))||(~isreal(parameters.pulse_dur))||...
       (numel(parameters.pulse_dur)~=1)||(parameters.pulse_dur<=0)
    error('parameters.pulse_dur should be a positive real scalar.');
end
if ~isfield(parameters,'pulse_amp')
    error('pulse amplitude must be supplied in parameters.pulse_amp field.');
elseif (~isnumeric(parameters.pulse_amp))||(~isreal(parameters.pulse_amp))||...
       (numel(parameters.pulse_amp)~=1)
    error('parameters.pulse_amp should be a real scalar.');
end
if ~isfield(parameters,'pulse_num')
    error('number of pulses must be supplied in parameters.pulse_num field.');
elseif (~isnumeric(parameters.pulse_num))||(~isreal(parameters.pulse_num))||...
       (numel(parameters.pulse_num)~=1)||(parameters.pulse_num<1)||...
       (mod(parameters.pulse_num,1)~=0)
    error('parameters.pulse_num should be a positive integer.');
end
if ~isfield(parameters,'n_periods')
    error('number of rotor periods must be supplied in parameters.n_periods field.');
elseif (~isnumeric(parameters.n_periods))||(~isreal(parameters.n_periods))||...
       (numel(parameters.n_periods)~=1)||(parameters.n_periods<1)||...
       (mod(parameters.n_periods,1)~=0)
    error('parameters.n_periods should be a positive integer.');
end
end

% Many respectable physicists said that they weren't going to stand
% for this - partly because it was a debasement of science, but most-
% ly because they didn't get invited to those sorts of parties.
%
% Douglas Adams, "The Hitchhiker's Guide to the Galaxy"

