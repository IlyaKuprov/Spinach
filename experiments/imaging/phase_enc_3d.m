% 3D MRI pulse sequence, with a slice selection stage followed by phase-
% encoded acquisition of the slice. Syntax:
%
%          fid=phase_enc_3d(spin_system,parameters,H,R,K,G,F)
%     
% This sequence must be called from the imaging() context, which would
% provide H,R,K,G, and F. 
%
% Parameters: 
%
%  parameters.ss_grad_amp   -  the amplitude of slice selection 
%                              gradient,T/m
%
%  parameters.pe_grad_amp   -  phase encoding gradient amplitude, T/m
% 
%  parameters.ro_grad_amp   -  readout gradient amplitude, T/m
%
%  parameters.ss_grad_dur   -  the duration of the slice selection 
%                              gradient, seconds 
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
%  fid - k-space representation of the image 
%
% ilya.kuprov@weizmann.ac.uk
% ahmed.allami@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=phase_enc_3d.m>

function fid=phase_enc_3d(spin_system,parameters,H,R,K,G,F)

% Check consistency
grumble(spin_system,parameters,H,R,K,G,F);

% Assemble the background
B=H+F+1i*R+1i*K;

% Make pulse operators and kron them up into Fokker-Planck space
Sx=operator(spin_system,'Lx','1H'); Sx=kron(speye(prod(parameters.npts)),Sx);
Sy=operator(spin_system,'Ly','1H'); Sy=kron(speye(prod(parameters.npts)),Sy);

% Slice selection operator 
L=B+parameters.ss_grad_amp*G{1};

% Apply the slice selection pulse
parameters.rho0=shaped_pulse_af(spin_system,L,Sx,Sy,parameters.rho0,...
                                parameters.rf_frq_list,parameters.rf_amp_list,...
                                parameters.rf_dur_list,parameters.rf_phi,2,'expv');
                            
% Rollback operator
L=B-parameters.ss_grad_amp*G{1};                             
      
% Run the rollback gradient
parameters.rho0=evolution(spin_system,L,[],parameters.rho0,...
                          sum(parameters.rf_dur_list)/2,1,'final');
                      
% Evolve for the echo time
parameters.rho0=evolution(spin_system,B,[],parameters.rho0,parameters.t_echo,1,'final');

% Apply 180-degree pulse
parameters.rho0=step(spin_system,Sy,parameters.rho0,pi);

% Evolve for the echo time
parameters.rho0=evolution(spin_system,B,[],parameters.rho0,parameters.t_echo,1,'final');
     
% Project out Hx+1i*Hy in every voxel
mri_slice=fpl2phan(parameters.rho0,state(spin_system,'L+','1H'),parameters.npts);

% Get sample dimension information
dims=zeros(1,6); dims([1 3 5])=-parameters.dims; dims([2 4 6])=+parameters.dims;

% Draw the slice in three dimensions
figure(); volplot(abs(mri_slice),dims); ktitle('after slice sel. and echo'); drawnow();

% Get phase encoding gradient range
pe_grad_amps=linspace(-parameters.pe_grad_amp,...
                       parameters.pe_grad_amp,...
                       parameters.image_size(1));

% Get k-space sampling parameters
nsteps=parameters.image_size(2)-1; step_length=parameters.ro_grad_dur/nsteps;

% Preallocate the image
fid=zeros(parameters.image_size,'like',1i);

% Precompute evolution generators
F_preroll=B-parameters.ro_grad_amp*G{3};
F_readout=B+parameters.ro_grad_amp*G{3};

% Loop over frequencies
parfor n=1:parameters.image_size(1) %#ok<*PFBNS>
    
    % Run the phase encoding gradient
    rho=evolution(spin_system,B+pe_grad_amps(n)*G{2},[],parameters.rho0,parameters.pe_grad_dur,1,'final');
       
    % Run the pre-roll gradient
    rho=evolution(spin_system,F_preroll,[],rho,parameters.ro_grad_dur/2,1,'final');
    
    % Detect under the readout gradient
    fid(n,:)=evolution(spin_system,F_readout,parameters.coil,rho,step_length,nsteps,'observable');   

end

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
if ~isfield(parameters,'ss_grad_amp')
    error('slice selection gradient amplitude must be specified in parameters.ss_grad_amp field.');
end
if (~isnumeric(parameters.ss_grad_amp))||(~isreal(parameters.ss_grad_amp))||...
   (~isscalar(parameters.ss_grad_amp))
    error('parameters.ss_grad_amp must be a real scalar.');
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

% Right, as the world goes, is only in question between equals 
% in power, while the strong do what they can and the weak suf-
% fer what they must.
%
% Thucydides

