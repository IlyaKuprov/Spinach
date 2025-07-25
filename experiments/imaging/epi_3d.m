% Diffusion weighted 3D echo planar imaging pulse sequence. Syntax: 
% 
%           fid=epi_3d(spin_system,parameters,H,R,K,G,F)
%
% This sequence must be called from the imaging() context, which
% would provide H,R,K,G, and F. Parameters:
%
%  parameters.ss_grad_amp   -  the amplitude of slice selection 
%                              gradient,T/m
%
%  parameters.ss_grad_dur   -  the duration of the slice selection 
%                              gradient, seconds 
%
%  parameters.image_size    -  number of points in each dimension of
%                              the resulting image
%
%  parameters.pe_grad_amp   -  phase encoding gradient amplitude, T/m
% 
%  parameters.ro_grad_amp   -  readout gradient amplitude, T/m
%
%  parameters.grad_dur      -  the duration of the gradients, seconds
%
%  parameters.t_echo        -  echo time after slice selection, seconds
%
%  parameters.diff_g_amp    -  [optional] a vector of diffusion gra-
%                              dient pair amplitudes in X,Y (T/m) to
%                              be active during the echo time
%
% Outputs:
%
%  fid - k-space representation of the image 
% 
% ilya.kuprov@weizmann.ac.il
% ahmed.allami@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=epi_3d.m>

function fid=epi_3d(spin_system,parameters,H,R,K,G,F)

% Check consistency
grumble(spin_system,parameters,H,R,K,G,F);

% Assemble the background
B=H+F+1i*R+1i*K;

% Make pulse operators
Hx=operator(spin_system,'Lx','1H');
Hy=operator(spin_system,'Ly','1H');
Hx=kron(speye(prod(parameters.npts)),Hx);
Hy=kron(speye(prod(parameters.npts)),Hy);

% Apply the slice selection pulse
rho=shaped_pulse_af(spin_system,B+parameters.ss_grad_amp*G{1},...
                    Hx,Hy,parameters.rho0,parameters.rf_frq_list,...
                    parameters.rf_amp_list,parameters.rf_dur_list,...
                    parameters.rf_phi,2,'expv');
                            
% Run the rollback gradient
rho=evolution(spin_system,B-parameters.ss_grad_amp*G{1},...
              [],rho,sum(parameters.rf_dur_list)/2,1,'final');
                      
% Optional diffusion gradient
if isfield(parameters,'diff_g_amp')

    % Evolve for the echo time in the presence of diffusion gradient
    rho=step(spin_system,B+parameters.diff_g_amp(1)*G{1}+...
                           parameters.diff_g_amp(2)*G{2}+...
                           parameters.diff_g_amp(3)*G{3},...
             rho,parameters.t_echo);

else

    % Simply evolve for the echo time
    rho=evolution(spin_system,B,[],rho,parameters.t_echo,1,'final');

end

% Apply an ideal 180-degree pulse
rho=step(spin_system,Hy,rho,pi);

% Optional diffusion gradient
if isfield(parameters,'diff_g_amp')

    % Evolve for the echo time in the presence of diffusion gradient
    rho=step(spin_system,B+parameters.diff_g_amp(1)*G{1}+...
                           parameters.diff_g_amp(2)*G{2}+...
                           parameters.diff_g_amp(3)*G{3},...
             rho,parameters.t_echo);

else

    % Simply evolve for the echo time
    rho=evolution(spin_system,B,[],rho,parameters.t_echo,1,'final');

end
     
% Project out Hx+1i*Hy in every voxel
mri_slice=fpl2phan(rho,state(spin_system,'L+','1H'),parameters.npts);

% Get sample dimension information
dims=zeros(1,6); 
dims([1 3 5])=-parameters.dims/2; 
dims([2 4 6])=+parameters.dims/2;

% Draw the slice in three dimensions
figure(); volplot(abs(mri_slice),dims); 
ktitle('after slice sel. and echo'); drawnow();

% Preroll the gradients
rho=evolution(spin_system,B-parameters.pe_grad_amp*G{2}...
                           -parameters.ro_grad_amp*G{3},[],rho,...
                            parameters.pe_grad_dur/2,1,'final');

% Precompute propagators
P_ro_p=propagator(spin_system,B+parameters.ro_grad_amp*G{3},...
                  parameters.ro_grad_dur/(parameters.image_size(2)-1));
P_ro_m=propagator(spin_system,B-parameters.ro_grad_amp*G{3},...
                  parameters.ro_grad_dur/(parameters.image_size(2)-1));
P_pe_p=propagator(spin_system,B+parameters.pe_grad_amp*G{2},...
                  parameters.pe_grad_dur/(parameters.image_size(1)-1));

% Move to GPU if necessary
if ismember('gpu',spin_system.sys.enable)

    % Upload relevant objects
    rho=gpuArray(rho); P_ro_p=gpuArray(P_ro_p);
    P_ro_m=gpuArray(P_ro_m); P_pe_p=gpuArray(P_pe_p);
    coil=gpuArray(parameters.coil);

    % Preallocate k-space image
    fid=1i*gpuArray.zeros(parameters.image_size);

    % Inform the user
    report(spin_system,'running the EPI loop on GPU...');

else

    % Preallocate k-space image
    fid=zeros(parameters.image_size,'like',1i);

    % Get the coil
    coil=parameters.coil;

    % Inform the user
    report(spin_system,'running the EPI loop on CPU...');

end

% Phase encoding loop
for n=1:parameters.image_size(1)

    % Inform the user
    report(spin_system,['phase encoding rung ' int2str(n) ...
                        '/' int2str(parameters.image_size(1))]);
    
    % Determine readout gradient sign
    ro_grad_sign=2*mod(n,2)-1;
    
    % Readout loop
    for k=1:parameters.image_size(2)
        
        % Detect under readout gradient
        if ro_grad_sign>0
            fid(n,k)=hdot(coil,rho);
            rho=P_ro_p*rho;
        else
            rho=P_ro_m*rho;
            fid(n,parameters.image_size(2)-k+1)=hdot(coil,rho);
        end
        
    end
    
    % Propagate under encoding gradient
    rho=P_pe_p*rho;
    
end

% Retrieve the fid
fid=gather(fid);

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
if ~isfield(parameters,'ss_grad_dur')
    error(' slice selection gradient duration must be specified in parameters.ss_grad_dur field.');
end
if (~isnumeric(parameters.ss_grad_dur))||(~isreal(parameters.ss_grad_dur))||...
   (~isscalar(parameters.ss_grad_dur))||(parameters.ss_grad_dur<=0)
    error('parameters.ss_grad_dur must be a positive real scalar.');
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

% The problem with Australians is not that so many of them 
% are descended from convicts, but that so many are descen-
% ded from prison officers.
%
% Clive James

