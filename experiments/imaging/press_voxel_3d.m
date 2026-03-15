% Voxel selection diagnostics function for 3D PRESS sequences. Returns
% the sample excitation profile. Syntax:
%
%       phan=press_voxel_3d(spin_system,parameters,H,R,K,G,F)
%
% This sequence must be called from the imaging() context, which
% would provide H, R, K, G, and F. 
%
% Parameters:
%
%    parameters.ss_grad_amp - the three amplitudes of slice selection 
%                             gradient, T/m
%
%    parameters.rf_frq_list - cell array of three vectors of RF frequ-
%                             encies at each pulse slice, Hz
%
%    parameters.rf_amp_list - cell array of three vectors of RF 
%                             amplitudes at each pulse slice, rad/s
%
%    parameters.rf_dur_list - cell array of three vectors of pulse 
%                             slice durations, in seconds
%
%    parameters.rf_phi      - cell array of three pulse phases at 
%                             time zero
%
%    parameters.max_rank    - cell array of three maximum rank in the 
%                             Fokker-Planck pulse operator (2 is 
%                             usually enough)
%
% Outputs:
%
%    phan - the excitation profile imprinted into a 3D phantom.
%
% Notes: add 'polyadic' to sys.enable, or this will crash your computer.
%
% ilya.kuprov@weizmann.ac.il
% a.j.allami@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=press_voxel_3d.m>

function phan=press_voxel_3d(spin_system,parameters,H,R,K,G,F)

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

% X slice selection pulse
rho=shaped_pulse_af(spin_system,L+parameters.ss_grad_amp(1)*G{1},Lx,Ly,...
                    rho,parameters.rf_frq_list{1},parameters.rf_amp_list{1},...
                    parameters.rf_dur_list{1},parameters.rf_phi{1},...
                    parameters.max_rank{1});
            
% X rephasing gradient
rho=evolution(spin_system,L-parameters.ss_grad_amp(1)*G{1},[],...
              rho,sum(parameters.rf_dur_list{1})/2,1,'final');
            
% Isolate single-quantum
rho=coherence(spin_system,rho,{{parameters.spins{1},1}});
            
% Y slice selection pulse (scaled down to 90 degrees)
rho=shaped_pulse_af(spin_system,L+parameters.ss_grad_amp(2)*G{2},Lx,Ly,...
                    rho,parameters.rf_frq_list{2},parameters.rf_amp_list{2},...
                    parameters.rf_dur_list{2}/2,parameters.rf_phi{2},...
                    parameters.max_rank{2});
                
% Y rephasing gradient
rho=evolution(spin_system,L-parameters.ss_grad_amp(2)*G{2},[],...
              rho,sum(parameters.rf_dur_list{2})/4,1,'final');
          
% Isolate zero-quantum
rho=coherence(spin_system,rho,{{parameters.spins{1},0}});
          
% Z slice selection pulse (scaled down to 90 degrees)
rho=shaped_pulse_af(spin_system,L+parameters.ss_grad_amp(3)*G{3},Lx,Ly,...
                    rho,parameters.rf_frq_list{3},parameters.rf_amp_list{3},...
                    parameters.rf_dur_list{3}/4,parameters.rf_phi{3},...
                    parameters.max_rank{3});
                
% Z rephasing gradient
rho=evolution(spin_system,L-parameters.ss_grad_amp(3)*G{3},[],...
              rho,sum(parameters.rf_dur_list{3})/4,1,'final');
          
% Isolate single-quantum
rho=coherence(spin_system,rho,{{parameters.spins{1},1}});
                
% Get the phantom
Lp=state(spin_system,'L+',parameters.spins{1});
phan=abs(fpl2phan(rho,Lp,parameters.npts));
            
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
if (~iscell(G))||(numel(G)<3)
    error('the G argument must be a cell array with at least three gradient operators.');
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
   (numel(parameters.npts)~=3)||any(parameters.npts<1)||...
   any(mod(parameters.npts,1)~=0)
    error('parameters.npts must be a vector of three positive integers.');
end
if ~isfield(parameters,'ss_grad_amp')
    error('gradient amplitude must be specified in parameters.ss_grad_amp field.');
end
if (~isnumeric(parameters.ss_grad_amp))||(~isreal(parameters.ss_grad_amp))||...
   (numel(parameters.ss_grad_amp)~=3)
    error('parameters.ss_grad_amp must be a vector of three real numbers.');
end
if ~isfield(parameters,'rf_frq_list')
    error('parameters.rf_frq_list field must be present.');
end
if (~iscell(parameters.rf_frq_list))||(numel(parameters.rf_frq_list)~=3)||...
   any(cellfun(@(x) (~isnumeric(x))||(~isreal(x))||(~isvector(x)),parameters.rf_frq_list))
    error('parameters.rf_frq_list must be a cell array of three real vectors.');
end
if ~isfield(parameters,'rf_amp_list')
    error('parameters.rf_amp_list field must be present.');
end
if (~iscell(parameters.rf_amp_list))||(numel(parameters.rf_amp_list)~=3)||...
   any(cellfun(@(x) (~isnumeric(x))||(~isreal(x))||(~isvector(x)),parameters.rf_amp_list))
    error('parameters.rf_amp_list must be a cell array of three real vectors.');
end
if ~isfield(parameters,'rf_dur_list')
    error('parameters.rf_dur_list field must be present.');
end
if (~iscell(parameters.rf_dur_list))||(numel(parameters.rf_dur_list)~=3)||...
   any(cellfun(@(x) (~isnumeric(x))||(~isreal(x))||(~isvector(x))||any(x<=0),parameters.rf_dur_list))
    error('parameters.rf_dur_list must be a cell array of three positive real vectors.');
end
if ~isfield(parameters,'rf_phi')
    error('parameters.rf_phi field must be present.');
end
if (~iscell(parameters.rf_phi))||(numel(parameters.rf_phi)~=3)||...
   any(cellfun(@(x) (~isnumeric(x))||(~isreal(x))||(~isscalar(x)),parameters.rf_phi))
    error('parameters.rf_phi must be a cell array of three real scalars.');
end
if ~isfield(parameters,'max_rank')
    error('parameters.max_rank field must be present.');
end
if (~iscell(parameters.max_rank))||(numel(parameters.max_rank)~=3)||...
   any(cellfun(@(x) (~isnumeric(x))||(~isreal(x))||(~isscalar(x))||(x<1)||(mod(x,1)~=0),parameters.max_rank))
    error('parameters.max_rank must be a cell array of three positive integers.');
end
for n=1:3
    if (numel(parameters.rf_frq_list{n})~=numel(parameters.rf_amp_list{n}))||...
       (numel(parameters.rf_amp_list{n})~=numel(parameters.rf_dur_list{n}))
        error('RF frequency, amplitude, and duration lists must match within each slice-selection pulse.');
    end
end
end

% One of the hardest parts of working with academics is that
% you have been trained, over decades, to argue your case.
% Academics are stunningly persuasive and eloquent. All I can
% realistically do is pitch you against each other and hope 
% that my side wins... 
% 
% Southampton Chemistry 
% Department secretary, to IK

