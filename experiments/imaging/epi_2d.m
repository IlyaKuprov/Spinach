% Diffusion weighted echo planar 2D imaging pulse sequence with 
% variable diffusion encoding direction. Syntax:
% 
%           mri=epi_2d(spin_system,parameters,H,R,K,G,F)
% 
% This sequence must be called from the imaging() context, which
% would provide H,R,K,G, and F. Parameters: 
% 
%   parameters.pe_grad_dur  - the duration of the phase encoding
%                             gradient (X), seconds
% 
%   parameters.ro_grad_dur  - the duration of the readout 
%                             gradient (Y), seconds
%
%   parameters.image_size   - number of points in each dimension
%                             of the resulting image
%
%   parameters.diff_g_amp   - [optional] a vector of diffusion 
%                             gradient amplitudes in X,Y (T/m)
% 
%   parameters.diff_g_dur   - [optional] the duration of the 
%                             diffusion gradient, seconds
% 
% Outputs:
%
%   mri - MRI image with square sinebell apodisation.
% 
% a.j.allami@soton.ac.uk
% i.ka.kuprov@weizmann.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=epi_2d.m>

function mri=epi_2d(spin_system,parameters,H,R,K,G,F)

% Check consistency
grumble(spin_system,parameters,H,R,K,G,F);

% Assemble the background
B=H+F+1i*R+1i*K;

% Make pulse operators
Lx=operator(spin_system,'Lx',parameters.spins{1});
Ly=operator(spin_system,'Ly',parameters.spins{1});
Lx=kron(speye(prod(parameters.npts)),Lx);
Ly=kron(speye(prod(parameters.npts)),Ly);

% Apply an ideal 90-degree pulse
rho=step(spin_system,Ly,parameters.rho0,pi/2);

% Apply diffusion gradient
if isfield(parameters,'diff_g_amp')
    parameters.rho0=step(spin_system,B+parameters.diff_g_amp(1)*G{1}+...
                                       parameters.diff_g_amp(2)*G{2},...
                         parameters.rho0,parameters.diff_g_dur);
end

% Apply an ideal 180-degree pulse
parameters.rho0=step(spin_system,Lx,parameters.rho0,pi);

% Apply diffusion gradient
if isfield(parameters,'diff_g_amp')
    parameters.rho0=step(spin_system,B+parameters.diff_g_amp(1)*G{1}+...
                                       parameters.diff_g_amp(2)*G{2},...
                         parameters.rho0,parameters.diff_g_dur);
end

% Preroll the gradients
rho=evolution(spin_system,B-parameters.pe_grad_amp*G{1}...
                           -parameters.ro_grad_amp*G{2},[],rho,...
                            parameters.pe_grad_dur/2,1,'final');

% Preallocate k-space image
fid=zeros(parameters.image_size,'like',1i);

% Precompute propagators
P_ro_p=propagator(spin_system,B+parameters.ro_grad_amp*G{2},...
                  parameters.ro_grad_dur/(parameters.image_size(2)-1));
P_ro_m=propagator(spin_system,B-parameters.ro_grad_amp*G{2},...
                  parameters.ro_grad_dur/(parameters.image_size(2)-1));
P_pe_p=propagator(spin_system,B+parameters.pe_grad_amp*G{1},...
                  parameters.pe_grad_dur/(parameters.image_size(1)-1));

% Phase encoding loop
for n=1:parameters.image_size(1)
    
    % Determine readout gradient sign
    ro_grad_sign=2*mod(n,2)-1;
    
    % Readout loop
    for k=1:parameters.image_size(2)
        
        % Detect under readout gradient
        if ro_grad_sign>0
            fid(n,k)=parameters.coil'*rho;
            rho=P_ro_p*rho;
        else
            rho=P_ro_m*rho;
            fid(n,parameters.image_size(2)-k+1)=parameters.coil'*rho;
        end
        
    end
    
    % Propagate under encoding gradient
    rho=P_pe_p*rho;
    
end

% Apodisation
fid=apodisation(spin_system,fid,{{'sqsin'},{'sqsin'}});
 
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
if ~isfield(parameters,'ro_grad_dur')
    error('the frequency encoding gradient duration must be specified in parameters.ro_grad_dur field.');
end
if (~isnumeric(parameters.ro_grad_dur))||(~isreal(parameters.ro_grad_dur))||...
   (~isscalar(parameters.ro_grad_dur))||(parameters.ro_grad_dur<=0)
    error('parameters.ro_grad_dur must be a positive real scalar.');
end
if ~isfield(parameters,'pe_grad_dur')
    error('the phase encoding gradient duration must be specified in parameters.grad_dur field.');
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

% The man is dangerous - he believes what he says.
%
% Count de Mirabeau, about Robespierre

