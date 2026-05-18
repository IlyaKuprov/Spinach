% Horne-Morris gradient-selected COSY pulse sequence. Syntax:
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
%                            either 'P', 'N', or 'P+N', defaults
%                            to 'P'
%
%    H  - Hamiltonian matrix, received from context function
%
%    R  - relaxation superoperator, received from context function
%
%    K  - kinetics superoperator, received from context function
%
% Outputs:
%
%    fid  - two-dimensional free induction decay, or a structure
%           with P-type fid.pos and N-type fid.neg fields in
%           'P+N' mode
%
% Note: the default P-type pathway uses opposite gradient signs
%       and is less sensitive to mixing pulse phase errors. The
%       N-type pathway uses equal gradient signs.
%
% Note: 'P+N' mode returns P-type and N-type components for
%       echo/anti-echo recombination in phase-sensitive processing.
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

% Initial state up to a constant multiplier
rho=state(spin_system,'Lz',parameters.spins{1},'cheap');

% Detection state up to a constant multiplier
coil=state(spin_system,'L+',parameters.spins{1},'cheap');

% Get the pulse operator
Lx=operator(spin_system,'Lx',parameters.spins{1});

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

% Select the requested gradient pathway
switch upper(parameters.pathway)

    case 'P'

        % Use opposite signs for P-type selection
        rho_stack=grad_sandw(spin_system,L,rho_stack,P,parameters.g_amp*[1 -1],...
                             parameters.s_len,parameters.g_dur*[1 1],[1 1]);

        % Run the recovery delay after the second gradient
        if parameters.g_stab_del>0
            rho_stack=evolution(spin_system,L,[],rho_stack,parameters.g_stab_del,1,'final');
        end

        % Run the F2 evolution
        fid=evolution(spin_system,L,coil,rho_stack,timestep,...
                      parameters.npoints(2)-1,'observable');

    case 'N'

        % Use equal signs for N-type selection
        rho_stack=grad_sandw(spin_system,L,rho_stack,P,parameters.g_amp*[1 1],...
                             parameters.s_len,parameters.g_dur*[1 1],[1 1]);

        % Run the recovery delay after the second gradient
        if parameters.g_stab_del>0
            rho_stack=evolution(spin_system,L,[],rho_stack,parameters.g_stab_del,1,'final');
        end

        % Run the F2 evolution
        fid=evolution(spin_system,L,coil,rho_stack,timestep,...
                      parameters.npoints(2)-1,'observable');

    case 'P+N'

        % Select the P-type pathway with opposite gradient signs
        rho_stack_p=grad_sandw(spin_system,L,rho_stack,P,parameters.g_amp*[1 -1],...
                               parameters.s_len,parameters.g_dur*[1 1],[1 1]);

        % Select the N-type pathway with equal gradient signs
        rho_stack_n=grad_sandw(spin_system,L,rho_stack,P,parameters.g_amp*[1 1],...
                               parameters.s_len,parameters.g_dur*[1 1],[1 1]);

        % Run the recovery delay after the second gradient
        if parameters.g_stab_del>0
            rho_stack_p=evolution(spin_system,L,[],rho_stack_p,parameters.g_stab_del,1,'final');
            rho_stack_n=evolution(spin_system,L,[],rho_stack_n,parameters.g_stab_del,1,'final');
        end

        % Run the F2 evolution for both pathway components
        fid.pos=evolution(spin_system,L,coil,rho_stack_p,timestep,...
                          parameters.npoints(2)-1,'observable');
        fid.neg=evolution(spin_system,L,coil,rho_stack_n,timestep,...
                          parameters.npoints(2)-1,'observable');

end

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
elseif (~isnumeric(parameters.sweep))||(~isreal(parameters.sweep))||...
       (~isfinite(parameters.sweep))||(parameters.sweep<=0)
    error('parameters.sweep must be a positive real scalar.');
end
if ~isfield(parameters,'spins')
    error('working spins should be specified in parameters.spins variable.');
elseif numel(parameters.spins)~=1
    error('parameters.spins cell array should have exactly one element.');
elseif (~iscell(parameters.spins))||(~ischar(parameters.spins{1}))
    error('parameters.spins must be a one-element cell array of character strings.');
elseif ~ismember(parameters.spins{1},spin_system.comp.isotopes)
    error('parameters.spins refers to an isotope that is not present in the system.');
end
if ~isfield(parameters,'npoints')
    error('number of points should be specified in parameters.npoints variable.');
elseif numel(parameters.npoints)~=2
    error('parameters.npoints array should have exactly two elements.');
elseif (~isnumeric(parameters.npoints))||(~isreal(parameters.npoints))||...
       any(parameters.npoints<1)||any(mod(parameters.npoints,1)~=0)
    error('parameters.npoints must contain two positive integers.');
end
if ~isfield(parameters,'angle')
    error('second pulse angle should be specified in parameters.angle variable.');
elseif numel(parameters.angle)~=1
    error('parameters.angle array should have exactly one element.');
elseif (~isnumeric(parameters.angle))||(~isreal(parameters.angle))||...
       (~isfinite(parameters.angle))
    error('parameters.angle must be a finite real scalar.');
end
if ~isfield(parameters,'g_amp')
    error('gradient amplitude should be specified in parameters.g_amp variable.');
elseif (~isnumeric(parameters.g_amp))||(~isreal(parameters.g_amp))||...
       (~isscalar(parameters.g_amp))||(~isfinite(parameters.g_amp))||...
       (parameters.g_amp<=0)
    error('parameters.g_amp should be a positive real scalar.');
end
if ~isfield(parameters,'g_dur')
    error('gradient duration should be specified in parameters.g_dur variable.');
elseif (~isnumeric(parameters.g_dur))||(~isreal(parameters.g_dur))||...
       (~isscalar(parameters.g_dur))||(~isfinite(parameters.g_dur))||...
       (parameters.g_dur<=0)
    error('parameters.g_dur should be a positive real scalar.');
end
if ~isfield(parameters,'g_stab_del')
    error('gradient stabilisation delay should be specified in parameters.g_stab_del variable.');
elseif (~isnumeric(parameters.g_stab_del))||(~isreal(parameters.g_stab_del))||...
       (~isscalar(parameters.g_stab_del))||(~isfinite(parameters.g_stab_del))||...
       (parameters.g_stab_del<0)
    error('parameters.g_stab_del should be a non-negative real scalar.');
end
if ~isfield(parameters,'s_len')
    error('sample length should be specified in parameters.s_len variable.');
elseif (~isnumeric(parameters.s_len))||(~isreal(parameters.s_len))||...
       (~isscalar(parameters.s_len))||(~isfinite(parameters.s_len))||...
       (parameters.s_len<=0)
    error('parameters.s_len should be a positive real scalar.');
end
if ~isfield(parameters,'pathway')
    error('coherence pathway should be specified in parameters.pathway variable.');
elseif (~ischar(parameters.pathway))||...
       (~ismember(upper(parameters.pathway),{'P','N','P+N'}))
    error('parameters.pathway should be ''P'', ''N'', or ''P+N''.');
end
end

% For economic and political thought to make useful progress, 
% it needs to be informed by evolutionary biology. This seems
% a very necessary exercise, since any attempt to understand
% morality, politics, economics or business without reference
% to evolutionary biology is ridiculous. As I explain to my
% children, ants are Marxist, dogs are Burkean conservatives
% and cats are libertarians. And, as I explain to our clients,
% a flower is a weed with an advertising budget.
%
% Rory Sutherland
