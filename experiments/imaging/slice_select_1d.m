% Slice selection diagnostics function. Executes a shaped pulse on
% the user-supplied 1D phantom and records a 1D image. Syntax:
%
%       fid=slice_select_1d(spin_system,parameters,H,R,K,G,F)
%
% This sequence must be called from the imaging() context, which
% would provide H,R,K,G, and F. 
%
% Parameters:
%
%    parameters.ss_grad_amp - the amplitude of slice selection 
%                             gradient,T/m
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
%    fid - k-space signal, run an FT to get the image
%
% ilya.kuprov@weizmann.ac.il
% a.j.allami@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=slice_select_1d.m>

function fid=slice_select_1d(spin_system,parameters,H,R,K,G,F)

% Check consistency
grumble(spin_system,parameters,H,R,K,G,F);

% Compose Liouvillian
L=H+F+1i*R+1i*K;

% Get pulse operators
Lp=operator(spin_system,'L+',parameters.spins{1});
Lx=kron(speye(prod(parameters.npts)),(Lp+Lp')/2);
Ly=kron(speye(prod(parameters.npts)),(Lp-Lp')/2i);

% Slice selection pulse
parameters.rho0=shaped_pulse_af(spin_system,L+parameters.ss_grad_amp*G{1},Lx,Ly,...
                parameters.rho0,parameters.rf_frq_list,parameters.rf_amp_list,...
                parameters.rf_dur_list,parameters.rf_phi,parameters.max_rank);
            
% Rephasing gradient
parameters.rho0=evolution(spin_system,L-parameters.ss_grad_amp*G{1},[],...
                          parameters.rho0,sum(parameters.rf_dur_list)/2,1,'final');
            
% Call to hard imaging
fid=basic_1d_hard(spin_system,parameters,H,R,K,G,F);
            
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
if ~isfield(parameters,'ss_grad_amp')
    error('gradient amplitude must be specified in parameters.ss_grad_amp field.');
end
if (~isnumeric(parameters.ss_grad_amp))||(~isreal(parameters.ss_grad_amp))||...
   (~isscalar(parameters.ss_grad_amp))
    error('parameters.ss_grad_amp must be a real scalar.');
end
end

% "Beta" is Latin for "still doesn't work".
%
% Internet folklore

