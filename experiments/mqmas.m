% Rotor-synchronous MQMAS pulse sequence. Syntax:
%
%            fid=mqmas(spin_system,parameters,H,R,K)
%
% This function should normally be called using singlerot.m context
% that would provide H, R, and K.
%
% Parameters:
%
%     parameters.pulse_dur - duration of each pulse, a two-
%                            element vector, seconds
%
%     parameters.pulse_amp - amplitude of each pulse, a two-
%                            element vector, rad/s
%
%     parameters.mq_order  - MQMAS coherence order
%
%     parameters.rho0      - initial condition, usually Lz
%
%     parameters.coil      - detection state, usually L+
%
%     + the parameters required by the singlerot.m 
%       context function that will provide H, R, and K
%
% Outputs:
%
%     fid - 2D amplitude mode free induction decay
%
% Notes:
%
%     parameters.sweep should be set equal to parameters.rate
%
% i.kuprov@soton.ac.uk
% m.carravetta@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=mqmas.m>

function fid=mqmas(spin_system,parameters,H,R,K)

% Check consistency
grumble(spin_system,parameters,H,R,K);

% Compose Liouvillian
L=H+1i*R+1i*K;

% Get pulse operators
Lp=operator(spin_system,'L+',parameters.spins{1});
Lx=(Lp+Lp')/2; Lx=kron(speye(parameters.spc_dim),Lx);

% Apply the decoupling
[L,parameters.rho0]=decouple(spin_system,L,parameters.rho0,parameters.decouple);

% Run the first pulse
rho=step(spin_system,L+parameters.pulse_amp(1)*Lx,...
         parameters.rho0,parameters.pulse_dur(1));

% Do the coherence selection
rho=coherence(spin_system,rho,{{parameters.spins{1},parameters.mq_order}});

% Run the indirect dimension evolution
rho_stack=evolution(spin_system,L,[],rho,1/parameters.rate,...
                    parameters.npoints(1)-1,'trajectory');
                
% Run the second pulse
rho_stack=step(spin_system,L+parameters.pulse_amp(2)*Lx,...
               rho_stack,parameters.pulse_dur(2));

% Do the coherence selection
rho_stack=coherence(spin_system,rho_stack,{{parameters.spins{1},+1}});

% Run the direct dimension evolution
fid=evolution(spin_system,L,parameters.coil,rho_stack,1/parameters.rate,...
                            parameters.npoints(2)-1,'observable');

end

% Consistency enforcement
function grumble(spin_system,parameters,H,R,K)
if (~isnumeric(H))||(~isnumeric(R))||(~isnumeric(K))||...
   (~ismatrix(H))||(~ismatrix(R))||(~ismatrix(K))
    error('H, R and K arguments must be matrices.');
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
if abs(parameters.sweep)~=abs(parameters.rate)
    error('parameters.sweep must be equal to parameters.rate (rotor-synchronous sequence).');
end
end

% Trust not the telly, trust the fridge.
%
% A Russian saying

