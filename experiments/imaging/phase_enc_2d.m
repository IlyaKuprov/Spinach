% 2D phase encoding imaging pulse sequence with optional diffusion
% weighting during the echo time. Syntax: 
% 
%          mri=phase_enc(spin_system,parameters,H,R,K,G,F)
%
% This sequence must be called from the imaging() context, which
% would provide H,R,K,G, and F. 
%
% Parameters: 
%
%  parameters.t_echo        -  echo time, seconds
%
%  parameters.diff_g_amp   -   [optional] a vector of diffusion gra-
%                              dient pair amplitudes in X,Y (T/m) to
%                              be active during the echo time
%
%  parameters.pe_grad_amp   -  phase encoding gradient amplitude, T/m
% 
%  parameters.ro_grad_amp   -  readout gradient amplitude, T/m
%
%  parameters.pe_grad_dur   -  the duration of the phase encoding 
%                              gradient, seconds 
%
%  parameters.ro_grad_dur   -  the duration of the readout gradient,
%                              seconds
%
%  parameters.image_size    -  number of points in each dimension of
%                              the resulting image
%
% Outputs:
%
%  mri - MRI image with square sinebell apodisation
%
% i.kuprov@soton.ac.uk
% a.j.allami@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=phase_enc_2d.m>

function mri=phase_enc_2d(spin_system,parameters,H,R,K,G,F)

% Check consistency
grumble(spin_system,parameters,H,R,K,G,F);

% Assemble the background
B=H+F+1i*R+1i*K;

% Make pulse operators
Hy=operator(spin_system,'Ly',parameters.spins{1});
Hy=kron(speye(prod(parameters.npts)),Hy);

% Apply ideal 90-degree pulse
parameters.rho0=step(spin_system,Hy,parameters.rho0,pi/2);

% Optional diffusion encoding
if isfield(parameters,'diff_g_amp')
    
    % Evolve for the echo time under a diffusion gradient
    parameters.rho0=step(spin_system,B+parameters.diff_g_amp(1)*G{1}+...
                                       parameters.diff_g_amp(2)*G{2},...
                                       parameters.rho0,parameters.t_echo);
else

    % Just evolve for the echo time
    parameters.rho0=step(spin_system,B,parameters.rho0,parameters.t_echo);

end

% Apply 180-degree pulse
parameters.rho0=step(spin_system,Hy,parameters.rho0,pi);

% Optional diffusion encoding
if isfield(parameters,'diff_g_amp')
    
    % Evolve for the echo time under a diffusion gradient
    parameters.rho0=step(spin_system,B+parameters.diff_g_amp(1)*G{1}+...
                                       parameters.diff_g_amp(2)*G{2},...
                                       parameters.rho0,parameters.t_echo);
else

    % Just evolve for the echo time
    parameters.rho0=step(spin_system,B,parameters.rho0,parameters.t_echo);

end

% Get PE gradient range
pe_grad_amps=linspace(-parameters.pe_grad_amp,...
                      +parameters.pe_grad_amp,...
                       parameters.image_size(1));

% Preallocate the image
fid=zeros(parameters.image_size,'like',1i);

% Phase encoding loop
parfor n=1:parameters.image_size(1) %#ok<*PFBNS>
    
    % Phase encoding gradient
    rho=evolution(spin_system,B+pe_grad_amps(n)*G{1},[],parameters.rho0,parameters.pe_grad_dur,1,'final');
       
    % Get the timing parameters
    nsteps=parameters.image_size(2)-1; step_length=parameters.ro_grad_dur/nsteps;
    
    % Run the pre-roll gradient
    rho=evolution(spin_system,B-parameters.ro_grad_amp*G{2},[],rho,step_length/2,nsteps,'final');
    
    % Detect under the readout gradient
    fid(n,:)=evolution(spin_system,B+parameters.ro_grad_amp*G{2},parameters.coil,rho,step_length,nsteps,'observable');   

end

% Apodization
fid=apodization(fid,'sqsinbell-2d');

% Fourier transform
mri=real(fftshift(fft2(ifftshift(fid))));

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
if ~isfield(parameters,'t_echo')
    error('echo time must be specified in parameters.t_echo field.');
end
if (~isnumeric(parameters.t_echo))||(~isreal(parameters.t_echo))||...
   (~isscalar(parameters.t_echo))||(parameters.t_echo<=0)
    error('parameters.t_echo must be a positive real scalar.');
end
if ~isfield(parameters,'ro_grad_dur')
    error(' readout gradient duration must be specified in parameters.ro_grad_dur field.');
end
if (~isnumeric(parameters.ro_grad_dur))||(~isreal(parameters.ro_grad_dur))||...
   (~isscalar(parameters.ro_grad_dur))||(parameters.ro_grad_dur<=0)
    error('parameters.ro_grad_dur must be a positive real scalar.');
end
if ~isfield(parameters,'pe_grad_dur')
    error('phase encoding gradient duration must be specified in parameters.pe_grad_dur field.');
end
if (~isnumeric(parameters.pe_grad_dur))||(~isreal(parameters.pe_grad_dur))||...
   (~isscalar(parameters.pe_grad_dur))||(parameters.pe_grad_dur<=0)
    error('parameters.pe_grad_dur must be a positive real scalar.');
end
if ~isfield(parameters,'ro_grad_amp')
    error('readout gradient amplitude must be specified in parameters.ro_grad_amp field.');
end
if (~isnumeric(parameters.ro_grad_amp))||(~isreal(parameters.ro_grad_amp))||...
   (~isscalar(parameters.ro_grad_amp))
    error('parameters.ro_grad_amp must be a real scalar.');
end
if ~isfield(parameters,'pe_grad_amp')
    error('phase encoding gradient amplitude must be specified in parameters.pe_grad_amp field.');
end
if (~isnumeric(parameters.pe_grad_amp))||(~isreal(parameters.pe_grad_amp))||...
   (~isscalar(parameters.pe_grad_amp))
    error('parameters.pe_grad_amp must be a real scalar.');
end
end

% Morally wrong-footing rivals is one point of ideology, and once 
% everyone agrees on something (slavery is wrong) it ceases to be
% a significant moral issue because it no longer shows local riv-
% als in a bad light. Many argue that there are more slaves in the
% world today than in the 19th century. Yet because one's politi-
% cal rivals cannot be delegitimised by being on the wrong side of
% slavery, few care to be active abolitionists any more, compared
% to being, say, speech police.
%
% John Tooby

