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
if (~iscell(G))||(numel(G)<2)
    error('the G argument must be a cell array with at least two gradient operators.');
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
   (numel(parameters.npts)~=2)||any(parameters.npts<1)||...
   any(mod(parameters.npts,1)~=0)
    error('parameters.npts must be a vector of two positive integers.');
end
if ~isfield(parameters,'rho0')
    error('parameters.rho0 field must be present.');
end
if ~isnumeric(parameters.rho0)
    error('parameters.rho0 must be numeric.');
end
if ~isfield(parameters,'ss_grad_amp')
    error('gradient amplitude must be specified in parameters.ss_grad_amp field.');
end
if (~isnumeric(parameters.ss_grad_amp))||(~isreal(parameters.ss_grad_amp))||...
   (numel(parameters.ss_grad_amp)~=2)
    error('parameters.ss_grad_amp must be a vector of two real numbers.');
end
if ~isfield(parameters,'rf_frq_list')
    error('parameters.rf_frq_list field must be present.');
end
if (~iscell(parameters.rf_frq_list))||(numel(parameters.rf_frq_list)~=2)||...
   any(cellfun(@(x) (~isnumeric(x))||(~isreal(x))||(~isvector(x)),parameters.rf_frq_list))
    error('parameters.rf_frq_list must be a cell array of two real vectors.');
end
if ~isfield(parameters,'rf_amp_list')
    error('parameters.rf_amp_list field must be present.');
end
if (~iscell(parameters.rf_amp_list))||(numel(parameters.rf_amp_list)~=2)||...
   any(cellfun(@(x) (~isnumeric(x))||(~isreal(x))||(~isvector(x)),parameters.rf_amp_list))
    error('parameters.rf_amp_list must be a cell array of two real vectors.');
end
if ~isfield(parameters,'rf_dur_list')
    error('parameters.rf_dur_list field must be present.');
end
if (~iscell(parameters.rf_dur_list))||(numel(parameters.rf_dur_list)~=2)||...
   any(cellfun(@(x) (~isnumeric(x))||(~isreal(x))||(~isvector(x))||any(x<=0),parameters.rf_dur_list))
    error('parameters.rf_dur_list must be a cell array of two positive real vectors.');
end
if ~isfield(parameters,'rf_phi')
    error('parameters.rf_phi field must be present.');
end
if (~iscell(parameters.rf_phi))||(numel(parameters.rf_phi)~=2)||...
   any(cellfun(@(x) (~isnumeric(x))||(~isreal(x))||(~isscalar(x)),parameters.rf_phi))
    error('parameters.rf_phi must be a cell array of two real scalars.');
end
if ~isfield(parameters,'max_rank')
    error('parameters.max_rank field must be present.');
end
if (~iscell(parameters.max_rank))||(numel(parameters.max_rank)~=2)||...
   any(cellfun(@(x) (~isnumeric(x))||(~isreal(x))||(~isscalar(x))||(x<1)||(mod(x,1)~=0),parameters.max_rank))
    error('parameters.max_rank must be a cell array of two positive integers.');
end
if ~isfield(parameters,'sp_grad_amp')
    error('crusher gradient amplitude must be specified in parameters.sp_grad_amp field.');
end
if (~isnumeric(parameters.sp_grad_amp))||(~isreal(parameters.sp_grad_amp))||(~isscalar(parameters.sp_grad_amp))
    error('parameters.sp_grad_amp must be a real scalar.');
end
if ~isfield(parameters,'sp_grad_dur')
    error('crusher gradient duration must be specified in parameters.sp_grad_dur field.');
end
if (~isnumeric(parameters.sp_grad_dur))||(~isreal(parameters.sp_grad_dur))||...
   (~isscalar(parameters.sp_grad_dur))||(parameters.sp_grad_dur<=0)
    error('parameters.sp_grad_dur must be a positive real scalar.');
end
for n=1:2
    if (numel(parameters.rf_frq_list{n})~=numel(parameters.rf_amp_list{n}))||...
       (numel(parameters.rf_amp_list{n})~=numel(parameters.rf_dur_list{n}))
        error('RF frequency, amplitude, and duration lists must match within each slice-selection pulse.');
    end
end
end

% A great deal more was hidden in the Dirac equation than the author
% had expected when he wrote it down in 1928. Dirac himself remarked 
% in one of his talks that his equation was more intelligent than its
% author. It should be added, however, that it was Dirac who found
% most of the additional insights.
%
% Victor Weisskopf

