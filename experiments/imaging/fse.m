% Fast spin echo (FSE) pulse sequence. Syntax: 
% 
%             mri=fse(spin_system,parameters,H,R,K,G,F)
%
% This sequence must be called from the imaging() context, which
% would provide H,R,K,G, and F. 
%
% Parameters: 
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
%  mri - MRI image with square sinebell apodisation.
%
% a.j.allami@soton.ac.uk
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=fse.m>

function mri=fse(spin_system,parameters,H,R,K,G,F)

% Check consistency
grumble(spin_system,parameters,H,R,K,G,F);

% Assemble the background
B=H+F+1i*R+1i*K;

% Make pulse operators
Lp=operator(spin_system,'L+',parameters.spins{1});
Ly=kron(speye(prod(parameters.npts)),(Lp-Lp')/2i);

% Apply 90-degree pulse
rho=step(spin_system,Ly,parameters.rho0,pi/2);

% Move to left edge of k-space
rho=step(spin_system,B+parameters.ro_grad_amp*G{2},...
                     rho,parameters.ro_grad_dur/2);

% Preallocate k-space image
fid=zeros(parameters.image_size,'like',1i);

% Phase encoding gradient amplitudes
pe_grad_amps=linspace(-parameters.pe_grad_amp,...
                       parameters.pe_grad_amp,...
                       parameters.image_size(1));

% Get readout timing                   
timestep=parameters.ro_grad_dur/(parameters.image_size(2)-1);

% Phase encoding loop
for n=1:parameters.image_size(1)
    
    % 180-degree pulse
    rho=step(spin_system,Ly,rho,pi);
    
    % Kick up the k-space
    rho=step(spin_system,B+pe_grad_amps(n)*G{1},...
                         rho,parameters.pe_grad_dur);

    % Get trajectory using Krylov propagation
    rho_stack=krylov(spin_system,B+parameters.ro_grad_amp*G{2},...
                     parameters.coil,rho,timestep,...
                     parameters.image_size(2)-1,'trajectory');
    rho=rho_stack(:,end); fid(n,:)=(parameters.coil'*rho_stack)';
    
    % Kick down the k-space
    rho=step(spin_system,B-pe_grad_amps(n)*G{1},...
             rho,parameters.pe_grad_dur);
    
end

% Apodisation
fid=apodisation(spin_system,fid,{{'sqsin'},{'sqsin'}});
 
% Fourier transform
mri=-real(fftshift(fft2(ifftshift(fid)),2));

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

% I spent a lot of money on booze, birds and fast 
% cars. The rest I just squandered. 
%
% George Best

