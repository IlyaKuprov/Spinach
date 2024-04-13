% Spin echo pulse sequence. Syntax:
%
%        fid=spin_echo(spin_system,parameters,H,R,K,G,F)
%
% This sequence must be called from the imaging() context, which
% would provide H,R,K,G, and F. Parameters:
%
%      parameters.g_amp      - the amplitude of gradient T/m 
%
%      parameters.g_step_dur - time step duration
%
%      parameters.g_n_steps  - number of time steps
%
% Outputs:
%
%      fid - the time domain echo signal
%
% i.kuprov@soton.ac.uk
% a.j.allami@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=spin_echo.m>

function fid=spin_echo(spin_system,parameters,H,R,K,G,F)

% Check consistency
grumble(spin_system,parameters,H,R,K,G,F);

% Assemble the Liouvillian
L=H+F+1i*R+1i*K;

% Make pulse operators
Hp=operator(spin_system,'L+',parameters.spins{1});
Hy=kron(speye(parameters.npts),(Hp-Hp')/2i);

% Hard 90-degree pulse
rho=step(spin_system,Hy,parameters.rho0,pi/2);

% Evolution under the X gradient
rho=evolution(spin_system,L+parameters.g_amp*G{1},[],rho,parameters.g_step_dur,...
                                                         parameters.g_n_steps,'final');

% Hard 180-degree pulse
rho=step(spin_system,Hy,rho,pi);

% Detection under the X gradient
fid=evolution(spin_system,L+parameters.g_amp*G{1},parameters.coil,rho,...
                                                  parameters.g_step_dur,...
                                                2*parameters.g_n_steps,'observable');

end

% Consistency enforcement
function grumble(spin_system,parameters,H,R,K,G,F)
if ~strcmp(spin_system.bas.formalism,'sphten-liouv')
    error('this function is only available in sphten-liouv formalism.');
end
if (~isnumeric(H))||(~isnumeric(R))||(~isnumeric(K))||...
   (~isnumeric(F))||(~ismatrix(H))||(~ismatrix(R))||...
   (~ismatrix(K))||(~ismatrix(F))
    error('H,R,K,F arguments must be matrices.');
end
if (~all(size(H)==size(R)))||...
   (~all(size(R)==size(K)))||...
   (~all(size(K)==size(F)))
    error('H,R,K,F matrices must have the same dimension.');
end
if ~iscell(G)
    error('the G argument must be a cell array.');
end
if ~isfield(parameters,'g_amp')
    error('gradient amplitude must be specified in parameters.g_amp field.');
end
if (~isnumeric(parameters.g_amp))||(~isreal(parameters.g_amp))||...
   (~isscalar(parameters.g_amp))
    error('parameters.g_amp must be a real scalar.');
end
if ~isfield(parameters,'g_step_dur')
    error('gradient step duration must be specified in parameters.g_step_dur field.');
end
if (~isnumeric(parameters.g_step_dur))||(~isreal(parameters.g_step_dur))||...
   (~isscalar(parameters.g_step_dur))||(parameters.g_step_dur<=0)
    error('parameters.g_step_dur must be a positive real scalar.');
end
if ~isfield(parameters,'g_n_steps')
    error('number of gradient steps must be specified in parameters.g_n_steps field.');
end
if (~isnumeric(parameters.g_n_steps))||(~isreal(parameters.g_n_steps))||...
   (~isscalar(parameters.g_n_steps))||(mod(parameters.g_n_steps,1)~=0)||...
   (parameters.g_n_steps<1)
    error('parameters.g_n_steps must be a positive real integer.');
end
end

% Wicked thoughts and worthless efforts gradually set their 
% mark on the face, especially the eyes.
%
% Arthur Schopenhauer

