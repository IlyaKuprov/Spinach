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
% ilya.kuprov@weizmann.ac.il
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
   (~ismatrix(H))||(~ismatrix(R))||(~ismatrix(K))||...
   (~isequal(size(H),size(R),size(K)))||...
   (size(H,1)~=size(H,2))
    error('H, R and K arguments must be square matrices of equal size.');
end
if ~isfield(parameters,'spins')
    error('active spin must be specified in parameters.spins variable.');
end
if (~iscell(parameters.spins))||(numel(parameters.spins)~=1)||...
   (~all(cellfun(@ischar,parameters.spins)))
    error('parameters.spins must be a cell array containing a single isotope string.');
end
if ~ismember(parameters.spins{1},spin_system.comp.isotopes)
    error('the isotope specified in parameters.spins is not present in the system.');
end
if ~isfield(parameters,'pulse_dur')
    error('pulse durations must be specified in parameters.pulse_dur variable.');
end
if (~isnumeric(parameters.pulse_dur))||(~isreal(parameters.pulse_dur))||...
   (numel(parameters.pulse_dur)~=2)||any(parameters.pulse_dur<0)
    error('parameters.pulse_dur must be a two-element vector of non-negative real numbers.');
end
if ~isfield(parameters,'pulse_amp')
    error('pulse amplitudes must be specified in parameters.pulse_amp variable.');
end
if (~isnumeric(parameters.pulse_amp))||(~isreal(parameters.pulse_amp))||...
   (numel(parameters.pulse_amp)~=2)
    error('parameters.pulse_amp must be a two-element vector of real numbers.');
end
if ~isfield(parameters,'mq_order')
    error('multiple quantum order must be specified in parameters.mq_order variable.');
end
if (~isnumeric(parameters.mq_order))||(~isreal(parameters.mq_order))||...
   (~isscalar(parameters.mq_order))||(mod(parameters.mq_order,1)~=0)
    error('parameters.mq_order must be an integer.');
end
if ~isfield(parameters,'spc_dim')
    error('spatial problem dimension must be specified in parameters.spc_dim variable.');
end
if (~isnumeric(parameters.spc_dim))||(~isreal(parameters.spc_dim))||...
   (~isscalar(parameters.spc_dim))||(parameters.spc_dim<1)||...
   (mod(parameters.spc_dim,1)~=0)
    error('parameters.spc_dim must be a positive integer.');
end
if ~isfield(parameters,'npoints')
    error('number of points must be specified in parameters.npoints variable.');
end
if (~isnumeric(parameters.npoints))||(~isreal(parameters.npoints))||...
   (numel(parameters.npoints)~=2)||any(parameters.npoints<1)||...
   any(mod(parameters.npoints,1)~=0)
    error('parameters.npoints must be a two-element vector of positive integers.');
end
if ~isfield(parameters,'rate')
    error('MAS rate must be specified in parameters.rate variable.');
end
if (~isnumeric(parameters.rate))||(~isreal(parameters.rate))||...
   (~isscalar(parameters.rate))||(parameters.rate==0)
    error('parameters.rate must be a non-zero real scalar.');
end
if ~isfield(parameters,'sweep')
    error('sweep width must be specified in parameters.sweep variable.');
end
if (~isnumeric(parameters.sweep))||(~isreal(parameters.sweep))||...
   (~isscalar(parameters.sweep))||(parameters.sweep==0)
    error('parameters.sweep must be a non-zero real scalar.');
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

