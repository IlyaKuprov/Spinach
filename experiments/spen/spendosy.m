% Ultrafast DOSY pulse sequence. Syntax:
%
%             fid=spendosy(spin_system,parameters,H,R,K,G,F)
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
% parameters.offset         offset 
%
% parameters.cond           bondary conditions 
%
% parameters.Ga             acquisition gradient in T/m
%
% parameters.pulsenpoints   number of points in the pulse shape
%
% parameters.smfactor       smoothing factor for the pulse
%
% parameters.Te             duration of the pulse
%
% parameters.Tau            duration extra dephasing gradient
%
% parameters.BW             bandwidth of the pulse
%
% parameters.Ge             encoding gradient in T/m
%
% parameters.chirptype      can be 'wurst' or 'smoothed'
% 
% parameters.td             diffusion delay 
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
% <https://spindynamics.org/wiki/index.php?title=spendosy.m>

function fid=spendosy(spin_system,parameters,H,R,K,G,F)

% Check consistency
grumble(spin_system,parameters);

% Compose Liouvillian
L=H+F+1i*R+1i*K;

% Get pulse operators
Lp=operator(spin_system,'L+',parameters.spins{1});
Lx=kron(speye(prod(parameters.npts)),(Lp+Lp')/2);
Ly=kron(speye(prod(parameters.npts)),(Lp-Lp')/2i);

% Get the chirp pulse
[Cx,Cy]=chirp_pulse(parameters.pulsenpoints,parameters.Te,parameters.BW,...
                       parameters.smfactor,parameters.chirptype);

% Get the gradient amplitude list
gradient_amplitudes=parameters.Ge*ones(parameters.pulsenpoints,1);

% Apply the first pulse
rho=step(spin_system,Lx,parameters.rho0,pi/2);

% Select "+1" coherence
rho=coherence(spin_system,rho,{{parameters.spins{1},+1}});

% Apply chirp pulse with a positive gradient
time_grid=parameters.Te*ones(1,parameters.pulsenpoints)/parameters.pulsenpoints;
rho=shaped_pulse_xy(spin_system,L,{Lx,Ly,G{1}},{Cx,Cy,+gradient_amplitudes},...
                    time_grid,rho,'expv-pwc');

% Select "-1" coherence
rho=coherence(spin_system,rho,{{parameters.spins{1},-1}});

% Evolve under a positive gradient
rho=step(spin_system,L+parameters.Ge*G{1},rho,parameters.Tau);

% Apply the second pulse
rho=step(spin_system,Lx,rho,pi/2);
  
% Select "0" coherence
rho=coherence(spin_system,rho,{{parameters.spins{1},0}});

% Run the diffusion time evolution
rho=step(spin_system,L,rho,parameters.td-parameters.Tau-parameters.Te);

% Apply the third pulse 
rho=step(spin_system,Lx,rho,pi/2);

% Select "-1" coherence
rho=coherence(spin_system,rho,{{parameters.spins{1},-1}});

% Apply chirp pulse with a positive gradient
rho=shaped_pulse_xy(spin_system,L,{Lx,Ly,G{1}},{Cx,Cy,+gradient_amplitudes},...
                    time_grid,rho,'expv-pwc');

% Select "+1" coherence
rho=coherence(spin_system,rho,{{parameters.spins{1},+1}});

% Evolve under a positive gradient
rho=step(spin_system,L+parameters.Ge*G{1},rho,parameters.Tau);

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
    rho_stack{m}=gather(rho); 
    rho=PL*rho;
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
if ~isfield(parameters,'smfactor')
    error('smoothing factor should be specified in parameters.smfactor variable.');
elseif numel(parameters.smfactor)~=1
    error('parameters.smfactor array should have exactly one elements.');
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
if ~isfield(parameters,'chirptype')
    error('chiprtype, smoothed or wurst, should be specified in parameters.chirptype variable.');
end

end

% Critique is a description of what the critic would have done
% if he were the author. If he had the author's abilities, the
% author's audience, the author's luck, esteem, and personal 
% history. But he has none of that, and so he's a critic.
%
% Sergei Shnurov

