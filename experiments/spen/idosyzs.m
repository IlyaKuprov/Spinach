% A simplified model sequence of the ZS iDOSY pulse sequence. Syntax:
%
%           inten=idosyzs(spin_system,parameters,H,R,K,G,F)
% 
% This sequence must be called from the imaging() context, which
% would provide H, R, K, G and F. Parameters:
%
%    parameters.rho0           - initial state
%
%    parameters.coil           - detection state
%
%    parameters.spins          - nuclei on which the sequence runs  
%
%    parameters.g_amp          - gradient amplitude for diffusion 
%                                encoding (T/m)
%
%    parameters.sel_g_amp      - gradient amplitude during the 
%                                selective pulse (T/m)
%
%    parameters.rf_phi         - phase of the inversion pulse (rad/s)  
%                                        
%    parameters.delta_big      - length of the diffusion delay (s)
%
%    parameters.delta_sml      - pulse width of the gradient for 
%                                diffusion encoding (s)    
%
%    parameters.rf_dur         - length of the selective 180 
%                                degree pulse (s)
%
%    parameters.filename       - a character string specifying the 
%                                pulse shape
%
%    parameters.pulse_npoints  - number of points in the soft pulse
%
%    parameters.diff           - diffusion constant (m^2/s)
%
% Outputs:
%
%    inten  - the absolute value of the first point in
%             the free induction decay; this number is 
%             proportional to the integral of the real
%             part of the correctly phased spectrum
%
% mariagrazia.concilio@sjtu.edu.cn 
% gareth.morris@manchester.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=idosyzs.m>

function inten=idosyzs(spin_system,parameters,H,R,K,G,F)

% Check consistency
grumble(spin_system,parameters,H,R,K,F,G);

% Compose Liouvillian
L=H+F+1i*R+1i*K;

% Get pulse operators
Lp=operator(spin_system,'L+',parameters.spins{1});
Ly=kron(speye(prod(parameters.npts)),(Lp-Lp')/2i);
Lx=kron(speye(prod(parameters.npts)),(Lp+Lp')/2);

% Diffusion time evolution delay
parameters.delta=parameters.delta_big-parameters.delta_sml-parameters.rf_dur;

% Shaped pulse setup
[amps,phases,~,~,scaling_factor]=read_wave(parameters.filename,parameters.pulse_npoints);
time_grid=parameters.rf_dur*ones(1,parameters.pulse_npoints)/parameters.pulse_npoints;

% Gradient amplitudes for the gradient during the selective pulse
gradient_amplitudes=parameters.sel_g_amp*ones(parameters.pulse_npoints,1);

% Calculate the maximum RF field strength and calibrate the amplitude of
% the soft pulse (this calibration is the same as in BRUKER TopSpin)
gamma_B1=parameters.rf_phi/parameters.rf_dur;  % rad/s
gamma_B1_max=gamma_B1/scaling_factor;
amps=gamma_B1_max*amps;
         
% Pulse transformation from (amplitude, phase) to (X,Y)
[Cx,Cy]=polar2cartesian(amps,phases);

% Apply the first pulse
rho=step(spin_system,Ly,parameters.rho0,pi/2);

% Select coherence -1
rho=coherence(spin_system,rho,{{parameters.spins{1},-1}});

% Evolve under the first positive gradient
rho=evolution(spin_system,L+parameters.g_amp*G{1},[],rho,parameters.delta_sml,1,'final');

% Run the diffusion time evolution
rho=evolution(spin_system,L,[],rho,parameters.delta,1,'final');

% Execution of the 180 degree shaped pulse
rho=shaped_pulse_xy(spin_system,L,{Lx,Ly,G{1}},{Cx,Cy,+gradient_amplitudes},...
                    time_grid,rho,'expv-pwc');

% Select coherence +1
rho=coherence(spin_system,rho,{{parameters.spins{1},+1}});
    
% Evolve under second positive gradient
rho=evolution(spin_system,L+parameters.g_amp*G{1},[],rho,parameters.delta_sml,1,'final');

% Evolve under delay to obtain complete refocus of the signal 
rho=step(spin_system,L,rho,parameters.delta);

% Intensity is the first point in the FID
inten=abs(parameters.coil'*rho);

% Inform the user        
report(spin_system,['the scaling factor is ' num2str(scaling_factor)...
       ', the gamma B1 max is ' num2str(gamma_B1_max./(2*pi)) ' Hz ']);

end
 
% Consistency enforcement
function grumble(spin_system,parameters,H,R,K,F,G)
 if ~ismember(spin_system.bas.formalism,{'sphten-liouv'})
    error('this function is only available for sphten-liouv formalism.');
end
if (~isnumeric(H))||(~isnumeric(R))||(~isnumeric(K))||...
   (~isnumeric(F))||(~ismatrix(H))||(~ismatrix(R))||(~ismatrix(K))||(~ismatrix(F))
    error('H, R, K and F arguments must be matrices.');
end
if (~all(size(H)==size(R)))||(~all(size(R)==size(K)))||(~all(size(K)==size(F)))
    error('H, R, K and F matrices must have the same dimension.');
end
if ~iscell(G)
    error('the G must be a 1x3 cell array.');
end
if ~ismember(spin_system.bas.formalism,{'sphten-liouv'})
    error('this function is only available for sphten-liouv formalism.');
end
if ~isfield(parameters,'dims')
   error('sample dimension should be specified in parameters.dims variable.');
elseif numel(parameters.dims)~=1
    error('parameters.dims array should have exactly one element.');
end
if ~isfield(parameters,'npts')
    error('number of spin packets should be specified in parameters.npts variable.');
elseif numel(parameters.npts)~=1
    error('parameters.npts array should have exactly one element.');
end
if ~isfield(parameters,'spins')
    error('working spins should be specified in parameters.spins variable.');
elseif numel(parameters.spins)~=1
    error('parameters.spins cell array should have exactly one element.');
end
if ~isfield(parameters,'diff')
    error('the diffusion coefficient should be specified in parameters.diff variable.');
elseif numel(parameters.diff)~=1
    error('parameters.diff array should have exactly one element.');
end
if ~isfield(parameters,'g_amp')
    error('the gradient amplitude should be specified in parameters.g_amp variable.');
elseif numel(parameters.g_amp)~=1
    error('parameters.g_amp array should have exactly one element.');
end
 if ~isfield(parameters,'sel_g_amp')
    error('the amplitude of the gradient during the selective 180 degree pulse should be specified in parameters.sel_g_amp variable.');
elseif numel(parameters.sel_g_amp)~=1
   error('parameters.sel_g_amp array should have exactly one element.'); 
end
if ~isfield(parameters,'rf_dur')
    error('the length of the selective 180 pulse should be specified in parameters.rf_dur variable.');
elseif numel(parameters.rf_dur)~=1
    error('parameters.rf_dur array should have exactly one element.');
end
if ~isfield(parameters,'rf_phi')
    error('the phase of the selective pulse has to be specified in arameters.rf_phi variable.');
elseif numel(parameters.rf_phi)~=1
    error('parameters.rf_phi array should have exactly one element.');
end
if ~isfield(parameters,'delta_big')
    error('the diffusion delay should be specified in parameters.delta_big variable.');
elseif numel(parameters.delta_big)~=1
    error('parameters.delta_big array should have exactly one element.'); 
end
if ~isfield(parameters,'delta_sml')
    error('the gradient pulse width should be specified in parameters.delta_sml variable.');
elseif numel(parameters.delta_sml)~=1
    error('parameters.delta_sml array should have exactly one element.');
end
end

 