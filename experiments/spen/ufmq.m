% Ultrafast multiple-quantum NMR, a literal implementation of Figure 1A
% from (http://dx.doi.org/10.1002/cphc.201800667). Syntax:
%
%             fid=ufmq_nmr(spin_system,parameters,H,R,K,G,F)
%
% This sequence must be called from the imaging() context, which
% would provide H, R, K, G, and F. Parameters:
%
% parameters.spins          nuclei on which the sequence runs
%
% parameters.dims           size of the sample, m
%
% parameters.npts           number of grid points 
%
% parameters.npoints        number of acquired points for each 
%                           gradient readout
%
% parameters.nloops         number of loop, where each loop consists of 
%                           a positive and a negative readout
%
% parameters.offset         offset, Hz 
%
% parameters.Ga             acquisition gradient in T/m
%
% parameters.deltat         timestep for acquisition, s
%
% parameters.pulsenpoints   number of points in the pulse shape
%
% parameters.BW             bandwidth of the pulse, Hz
%
% parameters.Ge             encoding gradient, T/m
%
% parameters.Te             duration of the pulse, s
%
% parameters.chirptype      can be 'wurst' or 'smoothed'
%
% parameters.smfactor       smoothing factor for the pulse
%
% Outputs:
%
%  fid -  free induction decay of the ultrafast NMR spectrum.
% 
% m.g.concilio@soton.ac.uk
% jean-nicolas.dumez@univ-nantes.fr
% ilya.kuprov@weizmann.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=ufmq.m>

function fid=ufmq(spin_system,parameters,H,R,K,G,F)

% Check consistency
grumble(spin_system,parameters);

% Compose Liouvillian
L=H+F+1i*R+1i*K;

% Get pulse operators
Lp=operator(spin_system,'L+',parameters.spins{1});
Lx=polyadic({{speye(prod(parameters.npts)),(Lp+Lp')/2}});
Ly=polyadic({{speye(prod(parameters.npts)),(Lp-Lp')/2i}});
if ~ismember('polyadic',spin_system.sys.enable)
    Lx=inflate(Lx); Ly=inflate(Ly);
end

% Get chirp pulse waveform
[Cx,Cy]=chirp_pulse(parameters.pulsenpoints,parameters.Te,...
                       parameters.BW,parameters.nWURST,parameters.chirptype);
time_grid=parameters.Te*ones(1,parameters.pulsenpoints)/parameters.pulsenpoints;

% Get the gradient pulse waveform for spatial encoding
gradient_amplitudes=parameters.Ge(1)*ones(parameters.pulsenpoints,1);

% Apply the first pulse
rho=step(spin_system,Lx,parameters.rho0,pi/2);

% Apply first delay
rho=evolution(spin_system,L,[],rho,parameters.delay,1,'final');

% Apply 180 deg pulse
rho=step(spin_system,Lx,rho,pi);

% Apply second delay
rho=evolution(spin_system,L,[],rho,parameters.delay,1,'final');

% Select operator for the second 90 deg pulse
if mod(parameters.mqorder,2)==0
    
    % Apply the second 90 deg pulse on x
    rho=step(spin_system,Lx,rho,pi/2);    
    
else
    
    % Apply the second 90 deg pulse on y
    rho=step(spin_system,Ly,rho,pi/2);        
    
end

% Apply the pair of chirp pulses with opposite gradients
rho=coherence(spin_system,rho,{{parameters.spins{1},+parameters.mqorder}});
rho=shaped_pulse_xy(spin_system,L,{Lx,Ly,G{1}},{Cx,Cy,+gradient_amplitudes},...
                    time_grid,rho,'expv-pwc');
rho=coherence(spin_system,rho,{{parameters.spins{1},-parameters.mqorder}});
rho=shaped_pulse_xy(spin_system,L,{Lx,Ly,G{1}},{Cx,Cy,-gradient_amplitudes},...
                    time_grid,rho,'expv-pwc');
rho=coherence(spin_system,rho,{{parameters.spins{1},+parameters.mqorder}});

% Apply third 90 deg pulse
rho=step(spin_system,Lx,rho,pi/2);

% Coherence selection
rho=coherence(spin_system,rho,{{parameters.spins{1},+1}});

% Apply prephasing (to set the echo in the right position)
Taq=parameters.npoints*parameters.deltat;
rho=step(spin_system,L-parameters.Ga*G{1},rho,Taq/2);

% Preallocate the fid
fid=zeros([parameters.npoints parameters.nloops],'like',1i);

% Upload to GPU
if ismember('gpu',spin_system.sys.enable)
    L=gpuArray(L); G{1}=gpuArray(G{1}); 
    rho=gpuArray(rho); fid=gpuArray(fid);
end

% Build evolution generators
EGp=L+parameters.Ga*G{1};
EGm=L-parameters.Ga*G{1};

% SPEN readout loop
for m=1:parameters.nloops
    
    % Tell the user
    report(spin_system,['readout loop ' num2str(m) '/' num2str(parameters.nloops) '...']);
    
    % Evolution under positive grad
    rho=step(spin_system,EGp,rho,Taq);
    
    % Detection under negative grad
    for k=1:parameters.npoints
        
        % Record magnetisation
        fid(k,m)=parameters.coil'*rho;
        
        % Make a time step
        rho=step(spin_system,EGm,rho,parameters.deltat);
        
    end
    
end

% Retrieve from GPU
if ismember('gpu',spin_system.sys.enable)
    fid=gather(fid);
end

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
if ~isfield(parameters,'mqorder')
    error('the multiple quantum coherence order should be specified in parameters.mqorder variable.');
end
end

% If many remedies are prescribed for an illness, you may
% be certain that the illness has no cure.
%
% Anton Chekhov

