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
%     parameters.rho0      - initial condition, usually Lz
%
%     parameters.coil      - detection state, usually L+
%
% Outputs:
%
%     fid                  - free induction decay
%
% ilya.kuprov@weizmann.ac.uk
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
if (~isnumeric(H))||(~isnumeric(R))||(~isnumeric(K))||(~ismatrix(H))||(~ismatrix(R))||(~ismatrix(K))
    error('H, R and K arguments must be matrices.');
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

