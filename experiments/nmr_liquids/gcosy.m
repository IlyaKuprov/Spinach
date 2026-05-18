% Gradient-selected COSY pulse sequence. Syntax:
%
%               fid=gcosy(spin_system,parameters,H,R,K)
%
% Parameters:
%
%    parameters.sweep        sweep width in Hz
%
%    parameters.npoints      number of points for both dimensions
%
%    parameters.spins        nuclei on which the sequence runs,
%                            specified as {'1H'}, {'13C'}, etc.
%
%    parameters.angle        second pulse angle in radians, usu-
%                            ally pi/2, but also allows COSY45,
%                            COSY60, etc.
%
%    parameters.g_amp        gradient amplitude in Gauss/cm,
%                            defaults to 3
%
%    parameters.g_dur        gradient duration in seconds,
%                            defaults to 2e-3
%
%    parameters.g_stab_del   post-gradient stabilisation delay in
%                            seconds, defaults to 2e-4
%
%    parameters.s_len        active sample length in cm,
%                            defaults to 1.5
%
%    parameters.pathway      optional coherence pathway selection,
%                            either 'P' or 'N', defaults to 'P'
%
%    H  - Hamiltonian matrix, received from context function
%
%    R  - relaxation superoperator, received from context function
%
%    K  - kinetics superoperator, received from context function
%
% Outputs:
%
%    fid  - two-dimensional free induction decay
%
% Note: the default P-type pathway uses opposite gradient signs
%       and is less sensitive to mixing pulse phase errors. The
%       N-type pathway uses equal gradient signs.
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=gcosy.m>

function fid=gcosy(spin_system,parameters,H,R,K)

% Set default gradient amplitude
if ~isfield(parameters,'g_amp')
    parameters.g_amp=3;
end

% Set default gradient duration
if ~isfield(parameters,'g_dur')
    parameters.g_dur=2e-3;
end

% Set default gradient stabilisation delay
if ~isfield(parameters,'g_stab_del')
    parameters.g_stab_del=2e-4;
end

% Set default active sample length
if ~isfield(parameters,'s_len')
    parameters.s_len=1.5;
end

% Set default coherence pathway
if ~isfield(parameters,'pathway')
    parameters.pathway='P';
end

% Check consistency
grumble(spin_system,parameters,H,R,K);

% Compose Liouvillian
L=H+1i*R+1i*K;

% Compute the evolution timestep
timestep=1/parameters.sweep;

% Initial state
rho=state(spin_system,'Lz',parameters.spins{1},'cheap');

% Detection state
coil=state(spin_system,'L+',parameters.spins{1},'cheap');

% Get the pulse operator
Lp=operator(spin_system,'L+',parameters.spins{1});
Lx=(Lp+Lp')/2;

% Apply the first pulse
rho=step(spin_system,Lx,rho,pi/2);

% Run the F1 evolution
rho_stack=evolution(spin_system,L,[],rho,timestep,parameters.npoints(1)-1,'trajectory');

% Build the second pulse propagator
P=propagator(spin_system,Lx,parameters.angle);

% Include the recovery delay after the first gradient
if parameters.g_stab_del>0
    P=P*propagator(spin_system,L,parameters.g_stab_del);
end

% Set the gradient signs
switch upper(parameters.pathway)

    case 'P'

        % Use opposite signs for P-type selection
        grad_signs=[1 -1];

    case 'N'

        % Use equal signs for N-type selection
        grad_signs=[1 1];

end

% Select the coherence pathway with a gradient sandwich
rho_stack=grad_sandw(spin_system,L,rho_stack,P,parameters.g_amp*grad_signs,...
                     parameters.s_len,parameters.g_dur*[1 1],[1 1]);

% Run the recovery delay after the second gradient
if parameters.g_stab_del>0
    rho_stack=evolution(spin_system,L,[],rho_stack,parameters.g_stab_del,1,'final');
end

% Run the F2 evolution
fid=evolution(spin_system,L,coil,rho_stack,timestep,parameters.npoints(2)-1,'observable');

end

% Consistency enforcement
function grumble(spin_system,parameters,H,R,K)
if ~ismember(spin_system.bas.formalism,{'sphten-liouv'})
    error('this function is only available for sphten-liouv formalism.');
end
if (~isnumeric(H))||(~isnumeric(R))||(~isnumeric(K))||...
   (~ismatrix(H))||(~ismatrix(R))||(~ismatrix(K))
    error('H, R and K arguments must be matrices.');
end
if (~all(size(H)==size(R)))||(~all(size(R)==size(K)))
    error('H, R and K matrices must have the same dimension.');
end
if ~isfield(parameters,'sweep')
    error('sweep width should be specified in parameters.sweep variable.');
elseif numel(parameters.sweep)~=1
    error('parameters.sweep array should have exactly one element.');
end
if ~isfield(parameters,'spins')
    error('working spins should be specified in parameters.spins variable.');
elseif numel(parameters.spins)~=1
    error('parameters.spins cell array should have exactly one element.');
end
if ~isfield(parameters,'npoints')
    error('number of points should be specified in parameters.npoints variable.');
elseif numel(parameters.npoints)~=2
    error('parameters.npoints array should have exactly two elements.');
end
if ~isfield(parameters,'angle')
    error('second pulse angle should be specified in parameters.angle variable.');
elseif numel(parameters.angle)~=1
    error('parameters.angle array should have exactly one element.');
end
if ~isfield(parameters,'g_amp')
    error('gradient amplitude should be specified in parameters.g_amp variable.');
elseif (~isnumeric(parameters.g_amp))||(~isreal(parameters.g_amp))||...
       (~isscalar(parameters.g_amp))||(parameters.g_amp<=0)
    error('parameters.g_amp should be a positive real scalar.');
end
if ~isfield(parameters,'g_dur')
    error('gradient duration should be specified in parameters.g_dur variable.');
elseif (~isnumeric(parameters.g_dur))||(~isreal(parameters.g_dur))||...
       (~isscalar(parameters.g_dur))||(parameters.g_dur<=0)
    error('parameters.g_dur should be a positive real scalar.');
end
if ~isfield(parameters,'g_stab_del')
    error('gradient stabilisation delay should be specified in parameters.g_stab_del variable.');
elseif (~isnumeric(parameters.g_stab_del))||(~isreal(parameters.g_stab_del))||...
       (~isscalar(parameters.g_stab_del))||(parameters.g_stab_del<0)
    error('parameters.g_stab_del should be a non-negative real scalar.');
end
if ~isfield(parameters,'s_len')
    error('sample length should be specified in parameters.s_len variable.');
elseif (~isnumeric(parameters.s_len))||(~isreal(parameters.s_len))||...
       (~isscalar(parameters.s_len))||(parameters.s_len<=0)
    error('parameters.s_len should be a positive real scalar.');
end
if ~isfield(parameters,'pathway')
    error('coherence pathway should be specified in parameters.pathway variable.');
elseif (~ischar(parameters.pathway))||(~ismember(upper(parameters.pathway),{'P','N'}))
    error('parameters.pathway should be either ''P'' or ''N''.');
end
end

% Bring me into the company of those who seek truth, and deliver me from
% those who have found it.
%
% Arthur C. Clarke

