% Voxel selection diagnostics function for 1D PRESS sequences. Re-
% turns the sample excitation profile. Syntax:
%
%       phan=press_voxel_1d(spin_system,parameters,H,R,K,G,F)
%
% This sequence must be called from the imaging() context, which
% would provide H,R,K,G, and F. 
%
% Parameters:
%
%    parameters.ss_grad_amp - the amplitude of slice selection 
%                             gradient, T/m
%
%    parameters.rf_frq_list - a vector of RF frequencies at each
%                             pulse slice, Hz
%
%    parameters.rf_amp_list - a vector of RF amplitudes at each
%                             pulse slice, rad/s
%
%    parameters.rf_dur_list - a vector of pulse slice durations,
%                             in seconds
%
%    parameters.rf_phi      - pulse phase at time zero
%
%    parameters.max_rank    - maximum rank in the Fokker-Planck
%                             pulse operator (2 is usually enough)
%
% Outputs:
%
%    phan - the excitation profile imprinted into a 1D phantom.
%
% ilya.kuprov@weizmann.ac.il
% a.j.allami@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=press_voxel_1d.m>

function phan=press_voxel_1d(spin_system,parameters,H,R,K,G,F)

% Check consistency
grumble(spin_system,parameters,H,R,K,G,F);

% Compose Liouvillian
L=H+F+1i*R+1i*K;

% Get pulse operators
Lp=operator(spin_system,'L+',parameters.spins{1});
Lx=polyadic({{opium(prod(parameters.npts),1),(Lp+Lp')/2}});
Ly=polyadic({{opium(prod(parameters.npts),1),(Lp-Lp')/2i}});

% Override input with uniform initial condition
Lz=state(spin_system,'Lz',parameters.spins{1});
rho=kron(ones(prod(parameters.npts),1),Lz);

% Slice selection pulse
rho=shaped_pulse_af(spin_system,L+parameters.ss_grad_amp*G{1},Lx,Ly,...
                    rho,parameters.rf_frq_list,parameters.rf_amp_list,...
                    parameters.rf_dur_list,parameters.rf_phi,...
                    parameters.max_rank);
            
% Rephasing gradient
rho=evolution(spin_system,L-parameters.ss_grad_amp*G{1},[],...
              rho,sum(parameters.rf_dur_list)/2,1,'final');
            
% Get the phantom
phan=fpl2phan(rho,Lz,[parameters.npts 1]);
            
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
if (~iscell(G))||(numel(G)<1)
    error('the G argument must be a cell array with at least one gradient operator.');
end
if ~isfield(parameters,'spins')
    error('parameters.spins field must be present.');
end
if (~iscell(parameters.spins))||(numel(parameters.spins)~=1)||(~ischar(parameters.spins{1}))
    error('parameters.spins must be a cell array containing one spin name.');
end
if ~isfield(parameters,'npts')
    error('parameters.npts field must be present.');
end
if (~isnumeric(parameters.npts))||(~isreal(parameters.npts))||...
   (~isscalar(parameters.npts))||(parameters.npts<1)||...
   (mod(parameters.npts,1)~=0)
    error('parameters.npts must be a positive integer.');
end
if ~isfield(parameters,'ss_grad_amp')
    error('gradient amplitude must be specified in parameters.ss_grad_amp field.');
end
if (~isnumeric(parameters.ss_grad_amp))||(~isreal(parameters.ss_grad_amp))||...
   (~isscalar(parameters.ss_grad_amp))
    error('parameters.ss_grad_amp must be a real scalar.');
end
if ~isfield(parameters,'rf_frq_list')
    error('parameters.rf_frq_list field must be present.');
end
if (~isnumeric(parameters.rf_frq_list))||(~isreal(parameters.rf_frq_list))||...
   (~isvector(parameters.rf_frq_list))
    error('parameters.rf_frq_list must be a real vector.');
end
if ~isfield(parameters,'rf_amp_list')
    error('parameters.rf_amp_list field must be present.');
end
if (~isnumeric(parameters.rf_amp_list))||(~isreal(parameters.rf_amp_list))||...
   (~isvector(parameters.rf_amp_list))
    error('parameters.rf_amp_list must be a real vector.');
end
if ~isfield(parameters,'rf_dur_list')
    error('parameters.rf_dur_list field must be present.');
end
if (~isnumeric(parameters.rf_dur_list))||(~isreal(parameters.rf_dur_list))||...
   (~isvector(parameters.rf_dur_list))||any(parameters.rf_dur_list<=0)
    error('parameters.rf_dur_list must be a positive real vector.');
end
if ~isfield(parameters,'rf_phi')
    error('parameters.rf_phi field must be present.');
end
if (~isnumeric(parameters.rf_phi))||(~isreal(parameters.rf_phi))||(~isscalar(parameters.rf_phi))
    error('parameters.rf_phi must be a real scalar.');
end
if ~isfield(parameters,'max_rank')
    error('parameters.max_rank field must be present.');
end
if (~isnumeric(parameters.max_rank))||(~isreal(parameters.max_rank))||...
   (~isscalar(parameters.max_rank))||(parameters.max_rank<1)||...
   (mod(parameters.max_rank,1)~=0)
    error('parameters.max_rank must be a positive integer.');
end
if (numel(parameters.rf_frq_list)~=numel(parameters.rf_amp_list))||...
   (numel(parameters.rf_amp_list)~=numel(parameters.rf_dur_list))
    error('parameters.rf_frq_list, parameters.rf_amp_list, and parameters.rf_dur_list must have the same number of elements.');
end
end

% Thank you for saving me from the Library. I can see 
% a lot of girls without boys there, and all they can
% do is study. It is a terrible fate.
% 
% One of IK's girlfriends

