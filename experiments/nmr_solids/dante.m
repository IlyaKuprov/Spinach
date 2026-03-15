% DANTE pulse sequence. Syntax:
%
%            fid=dante(spin_system,parameters,H,R,K)
%
% This function should normally be called using singlerot.m context
% that would provide H, R, and K.
%
% Parameters:
%
%     parameters.pulse_dur - duration of each pulse, seconds
%
%     parameters.pulse_amp - amplitude of each pulse, Hz
%
%     parameters.pulse_num - number of pulses within rotor period
%
%     parameters.n_periods - number of rotor periods that the 
%                            sequence is active for
%
%     parameters.spins     - working spin, specified as a
%                            single-element cell array
%
%     parameters.decouple  - isotopes to decouple, specified
%                            as a cell array
%
%     parameters.rate      - rotor frequency in Hz
%
%     parameters.sweep     - acquisition sweep width in Hz
%
%     parameters.npoints   - number of acquisition points
%
%     parameters.spc_dim   - Fokker-Planck spatial dimension
%
%     parameters.rho0      - initial condition, usually Lz
%
%     parameters.coil      - detection state, usually L+
%
% Outputs:
%
%     fid                  - free induction decay
%
% ilya.kuprov@weizmann.ac.il
% m.carravetta@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=dante.m>

function fid=dante(spin_system,parameters,H,R,K)

% Check consistency
grumble(spin_system,parameters,H,R,K);

% Compose Liouvillian
L=H+1i*R+1i*K;

% Get pulse operators
Lp=operator(spin_system,'L+',parameters.spins{1});
Lx=(Lp+Lp')/2; Lx=kron(speye(parameters.spc_dim),Lx);

% Apply the decoupling
[L,parameters.rho0]=decouple(spin_system,L,parameters.rho0,parameters.decouple);

% Timing parameters
rotor_period=abs(1/parameters.rate);
cycle_length=rotor_period/parameters.pulse_num;

% Bomb out if the schedule makes no sense
if (cycle_length-parameters.pulse_dur)<0
    error('DANTE pulse schedule does not fit into the rotor period.');
end

% Get the state vector going
rho=parameters.rho0;

% Loop over rotor cycles
for k=1:parameters.n_periods
    
    % Loop over pulses within the rotor cycle
    for n=1:parameters.pulse_num
        
        % Pulse
        rho=evolution(spin_system,L+2*pi*parameters.pulse_amp*Lx,...
                      [],rho,parameters.pulse_dur,1,'final');
        
        % Free evolution
        rho=evolution(spin_system,L,...
                      [],rho,cycle_length-parameters.pulse_dur,1,'final');
        
    end
    
end
    
% Run the evolution and watch the coil state
fid=evolution(spin_system,L,parameters.coil,rho,...
              1/parameters.sweep,parameters.npoints-1,'observable');

end

% Consistency enforcement
function grumble(spin_system,parameters,H,R,K)
if (~isnumeric(H))||(~isnumeric(R))||(~isnumeric(K))||...
   (~ismatrix(H))||(~ismatrix(R))||(~ismatrix(K))
    error('H, R and K arguments must be matrices.');
end
if (~all(size(H)==size(R)))||(~all(size(R)==size(K)))
    error('H, R and K matrices must have the same dimension.');
end
if ~isfield(parameters,'spins')
    error('working spin must be specified in parameters.spins variable.');
elseif (~iscell(parameters.spins))||(numel(parameters.spins)~=1)
    error('parameters.spins must be a single-element cell array.');
elseif ~all(ismember(parameters.spins,spin_system.comp.isotopes))
    error('parameters.spins must contain isotopes present in the system.');
end
if ~isfield(parameters,'spc_dim')
    error('Fokker-Planck dimension should be specified in parameters.spc_dim variable.');
elseif (~isnumeric(parameters.spc_dim))||(~isreal(parameters.spc_dim))||...
       (numel(parameters.spc_dim)~=1)||(parameters.spc_dim<1)||...
       (mod(parameters.spc_dim,1)~=0)
    error('parameters.spc_dim should be a positive integer.');
end
if ~isfield(parameters,'rate')
    error('MAS rate should be specified in parameters.rate variable.');
elseif (~isnumeric(parameters.rate))||(~isreal(parameters.rate))||...
       (numel(parameters.rate)~=1)||(parameters.rate==0)
    error('parameters.rate should be a non-zero real number.');
end
if ~isfield(parameters,'pulse_dur')
    error('pulse duration must be supplied in parameters.pulse_dur field.');
elseif (~isnumeric(parameters.pulse_dur))||(~isreal(parameters.pulse_dur))||...
       (numel(parameters.pulse_dur)~=1)||(parameters.pulse_dur<=0)
    error('parameters.pulse_dur should be a positive real number.');
end
if ~isfield(parameters,'pulse_amp')
    error('pulse amplitude must be supplied in parameters.pulse_amp field.');
elseif (~isnumeric(parameters.pulse_amp))||(~isreal(parameters.pulse_amp))||...
       (numel(parameters.pulse_amp)~=1)
    error('parameters.pulse_amp should be a real number.');
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
if ~isfield(parameters,'sweep')
    error('sweep width should be specified in parameters.sweep variable.');
end
if (~isnumeric(parameters.sweep))||(numel(parameters.sweep)~=1)||(~isreal(parameters.sweep))||(parameters.sweep<=0)
    error('parameters.sweep should be a positive real number.');
end
if ~isfield(parameters,'npoints')
    error('number of points should be specified in parameters.npoints variable.');
end
if (~isnumeric(parameters.npoints))||(numel(parameters.npoints)~=1)||...
   (~isreal(parameters.npoints))||(parameters.npoints<1)||(mod(parameters.npoints,1)~=0)
    error('parameters.npoints should be a positive integer.');
end
if ~isfield(parameters,'rho0')
    error('initial state must be specified in parameters.rho0 variable.');
end
if ~isfield(parameters,'coil')
    error('detection state must be specified in parameters.coil variable.');
end
if ~isfield(parameters,'decouple')
    error('list of decoupled spins (or an empty cell array) should be supplied in parameters.decouple variable.');
end
if ~iscell(parameters.decouple)
    error('parameters.decouple must be a cell array of strings.');
end
if numel(parameters.decouple)>0
    if any(~cellfun(@ischar,parameters.decouple))
        error('elements of parameters.decouple cell array must be strings.');
    end 
    if any(~ismember(parameters.decouple,spin_system.comp.isotopes))
        error('parameters.decouple contains isotopes that are not present in the system.');
    end
    if ~ismember(spin_system.bas.formalism,{'sphten-liouv'})
        error('analytical decoupling is only available for sphten-liouv formalism.');
    end
end
end

% Normality is a paved road: it is comfortable to walk,
% but no flowers grow on it.
%
% Vincent van Gogh

