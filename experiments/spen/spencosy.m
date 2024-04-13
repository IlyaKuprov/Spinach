% Ultrafast COSY pulse sequence. Syntax:
%
%         fid=spencosy(spin_system,parameters,H,R,K,G,F)
%
% The following parameters are required
%
% parameters.dims           size of the sample in m
%
% parameters.npts           number of spin packets
%
% parameters.spins          nuclei on which the sequence runs
%
% parameters.deltat         timestep for acquisition
%
% parameters.npoints        number of acquired points for each 
%                           gradient readout
%
% parameters.nloops         number of loop, where each loop consists of 
%                           a positive and a negative readout
%
% parameters.Ga             acquisition gradient in T/m
%
% parameters.pulsenpoints   number of points in the pulse shape
%
% parameters.nWURST         smoothing factor for the pulse
%
% parameters.Te             duration of the pulse
%
% parameters.BW             bandwidth of the pulse
%
% parameters.Ge             encoding gradient in T/m
%
% parameters.Gp             coherence selection gradient in T/m
%
% parameters.Tp             duration of the coherence selection gradient
%
% parameters.D              diffusion constant, m^2/s
%
% H                         Fokker-Planck Hamiltonian
%
% R                         Fokker-Planck relaxation superoperator
%
% K                         Fokker-Planck kinetics superoperator
%
% G                         Fokker-Planck gradient superoperators
%
% F                         Fokker-Planck diffusion and flow superoperator
%
% Note: the last five parameters are built automatically by the imaging
%       context function.
%
% jeannicolas.dumez@cnrs.fr
% i.kuprov@soton.ac.uk
% ludmilla.guduff@cnrs.fr
%
% <https://spindynamics.org/wiki/index.php?title=spencosy.m>

function fid=spencosy(spin_system,parameters,H,R,K,G,F)

% Check consistency
grumble(spin_system,parameters);

% Compose Liouvillian
L=H+F+1i*R+1i*K;

% Get pulse operators
Lp=operator(spin_system,'L+',parameters.spins{1});
Lx=kron(speye(prod(parameters.npts)),(Lp+Lp')/2);
Ly=kron(speye(prod(parameters.npts)),(Lp-Lp')/2i);

% Get chirp pulse waveform
[Cx,Cy]=chirp_pulse(parameters.pulsenpoints,parameters.Te,...
                       parameters.BW,parameters.nWURST,'wurst');

% Get the gradient pulse waveform
gradient_amplitudes=parameters.Ge*ones(parameters.pulsenpoints,1);

% Apply the first pulse
rho=step(spin_system,Lx,parameters.rho0,pi/2);

% Apply the pair of chirp pulses with opposite gradients
report(spin_system,'running chirps...')
time_grid=parameters.Te*ones(1,parameters.pulsenpoints)/parameters.pulsenpoints;
rho=shaped_pulse_xy(spin_system,L,{Lx,Ly,G{1}},{Cx,Cy,+gradient_amplitudes},...
                    time_grid,rho,'expv-pwc');
rho=shaped_pulse_xy(spin_system,L,{Lx,Ly,G{1}},{Cx,Cy,-gradient_amplitudes},...
                    time_grid,rho,'expv-pwc');

% Apply the second pulse with coherence selection
rho=step(spin_system,parameters.Gp*G{1},rho,parameters.Tp);
rho=step(spin_system,Lx,rho,pi/2);
rho=step(spin_system,parameters.Gp*G{1},rho,parameters.Tp);

% Apply prephasing
Taq=parameters.npoints*parameters.deltat;
rho=step(spin_system,L-parameters.Ga*G{1},rho,Taq/2);

% Build whole-loop propagators
P1L=propagator(spin_system,L+parameters.Ga*G{1},parameters.npoints*parameters.deltat);
P2L=propagator(spin_system,L-parameters.Ga*G{1},parameters.npoints*parameters.deltat);
PL=clean_up(spin_system,P2L*P1L,spin_system.tols.prop_chop); clear('P1L','P2L');

% Build the intra-loop propagator
P=propagator(spin_system,L+parameters.Ga*G{1},parameters.deltat);

% Move to the GPU if necessary
if ismember('gpu',spin_system.sys.enable)
    P=gpuArray(P); PL=gpuArray(PL);
    rho=gpuArray(full(rho)); 
    coil=gpuArray(parameters.coil);
else
    rho=full(rho); coil=parameters.coil;
end

% Generate loop starts
report(spin_system,'computing loop starts...');
rho_stack=cell(1,parameters.nloops); 
for m=1:parameters.nloops
    rho_stack{m}=rho; rho=PL*rho;
end

% Preallocate the fid
fid=zeros([parameters.npoints parameters.nloops],'like',1i);

% Propagate the system
report(spin_system,'computing loop bodies...');
parfor m=1:parameters.nloops %#ok<*PFBNS>
    rho=rho_stack{m};
    local_fid=zeros(parameters.npoints,1);
    for n=1:parameters.npoints 
        local_fid(n)=gather(coil'*rho); rho=P*rho;
    end
    fid(:,m)=local_fid;
end
report(spin_system,'propagation finished.');
          
end

% Consistency enforcement
function grumble(spin_system,parameters)
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
if ~isfield(parameters,'npoints')
    error('number of points should be specified in parameters.npoints variable.');
elseif numel(parameters.npoints)~=1
    error('parameters.npoints array should have exactly one elements.');
end
if ~isfield(parameters,'deltat')
    error('timestep should be specified in parameters.deltat variable.');
elseif numel(parameters.deltat)~=1
    error('parameters.deltat array should have exactly one elements.');
end
if ~isfield(parameters,'nloops')
    error('number of loops should be specified in parameters.nloops variable.');
elseif numel(parameters.nloops)~=1
    error('parameters.nloops array should have exactly one elements.');
end
if ~isfield(parameters,'Ga')
    error('acquisition gradient should be specified in parameters.Ga variable.');
elseif numel(parameters.Ga)~=1
    error('parameters.Ga array should have exactly one elements.');
end
if ~isfield(parameters,'pulsenpoints')
    error('number of points in the pulse shape should be specified in parameters.pulsenpoints variable.');
elseif numel(parameters.pulsenpoints)~=1
    error('parameters.pulsenpoints array should have exactly one elements.');
end
if ~isfield(parameters,'nWURST')
    error('smoothing factor should be specified in parameters.nWURST variable.');
elseif numel(parameters.nWURST)~=1
    error('parameters.nWURST array should have exactly one elements.');
end
if ~isfield(parameters,'Te')
    error('pulse duration should be specified in parameters.Te variable.');
elseif numel(parameters.Te)~=1
    error('parameters.Te array should have exactly one elements.');
end
if ~isfield(parameters,'BW')
    error('pulse bandwidth should be specified in parameters.BW variable.');
elseif numel(parameters.BW)~=1
    error('parameters.BW array should have exactly one elements.');
end
if ~isfield(parameters,'Ge')
    error('encoding gradient should be specified in parameters.Ge variable.');
elseif numel(parameters.Ge)~=1
    error('parameters.Ge array should have exactly one elements.');
end
if ~isfield(parameters,'Gp')
    error('coherence selection gradient should be specified in parameters.Gp variable.');
elseif numel(parameters.Gp)~=1
    error('parameters.Gp array should have exactly one elements.');
end
if ~isfield(parameters,'Tp')
    error('duration of coherence selection gradient should be specified in parameters.Tp variable.');
elseif numel(parameters.Tp)~=1
    error('parameters.Tp array should have exactly one elements.');
end
end

% Microsoft made a big deal about Windows NT getting a C2 security
% rating. They were much less forthcoming with the fact that this 
% rating only applied if the computer was not attached to a network
% and had no network card, and had its floppy drive epoxied shut,
% and was running on a Compaq 386. Solaris's C2 rating was just as
% silly.
%
% Bruce Schneier

