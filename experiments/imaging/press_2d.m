% 2D PRESS (voxel selective NMR) pulse sequence. Syntax:
%
%        fid=press_2d(spin_system,parameters,H,R,K,G,F)
%
% This sequence must be called from the imaging() context, which
% would provide H,R,K,G, and F. 
%
% Parameters:
%
%    parameters.ss_grad_amp - the two amplitudes of slice selection 
%                             gradient, T/m
%
%    parameters.rf_frq_list - cell array of two vectors of RF frequ-
%                             encies at each pulse slice, Hz
%
%    parameters.rf_amp_list - cell array of two vectors of RF 
%                             amplitudes at each pulse slice, rad/s
%
%    parameters.rf_dur_list - cell array of two vectors of pulse 
%                             slice durations, in seconds
%
%    parameters.rf_phi      - cell array of two pulse phases at 
%                             time zero
%
%    parameters.max_rank    - cell array of two maximum rank in the 
%                             Fokker-Planck pulse operator (2 is 
%                             usually enough)
%
%    parameters.sp_grad_amp - crusher gradient amplitude, T/m
%
%    parameters.sp_grad_dur - crusher gradient duration, seconds
%
% Outputs:
%
%    fid - free induction decay of the NMR spectrum
%
% ilya.kuprov@weizmann.ac.il
% a.j.allami@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=press_2d.m>

function fid=press_2d(spin_system,parameters,H,R,K,G,F)

% Check consistency
grumble(spin_system,parameters,H,R,K,G,F);

% Assemble the Liouvillian
L=H+F+1i*R+1i*K;

% Get pulse operators
Lp=operator(spin_system,'L+',parameters.spins{1});
Lx=kron(speye(prod(parameters.npts)),(Lp+Lp')/2);
Ly=kron(speye(prod(parameters.npts)),(Lp-Lp')/2i);

% Slice selection 90 degree pulse
parameters.rho0=shaped_pulse_af(spin_system,L+parameters.ss_grad_amp(1)*G{1},Lx,Ly,...
                parameters.rho0,parameters.rf_frq_list{1},parameters.rf_amp_list{1},...
                parameters.rf_dur_list{1},parameters.rf_phi{1},parameters.max_rank{1});
            
% Rephasing gradient
parameters.rho0=evolution(spin_system,L-parameters.ss_grad_amp(1)*G{1},[],...
                          parameters.rho0,sum(parameters.rf_dur_list{1})/2,1,'final');

% Apply a crusher gradient                       
parameters.rho0=evolution(spin_system,L+parameters.sp_grad_amp*(G{1}+G{2}),[],...
                          parameters.rho0,parameters.sp_grad_dur,1,'final');                      
                      
% Slice selection 180 degree pulse
parameters.rho0=shaped_pulse_af(spin_system,L+parameters.ss_grad_amp(2)*G{2},Lx,Ly,...
                parameters.rho0,parameters.rf_frq_list{2},2*parameters.rf_amp_list{2},...
                parameters.rf_dur_list{2},parameters.rf_phi{2},parameters.max_rank{2});
            
% Apply a crusher gradient                       
parameters.rho0=evolution(spin_system,L+parameters.sp_grad_amp*(G{1}+G{2}),[],...
                          parameters.rho0,parameters.sp_grad_dur,1,'final');   
               
% Acquisition
fid=acquire(spin_system,parameters,H+F,R,K);

end

% Consistency enforcement
function grumble(spin_system,parameters,H,R,K,G,F) %#ok<INUSL>
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
end

% A great deal more was hidden in the Dirac equation than the author
% had expected when he wrote it down in 1928. Dirac himself remarked 
% in one of his talks that his equation was more intelligent than its
% author. It should be added, however, that it was Dirac who found
% most of the additional insights.
%
% Victor Weisskopf

