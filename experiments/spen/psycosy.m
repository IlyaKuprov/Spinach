% Alan Kenwright's spatially encoded COSY sequence described in
% https://doi.org/10.1002/mrc.4727 (Figure 4). Syntax:
%
%         fid=psycosy(spin_system,parameters,H,R,K,G,F)
%
% Parameters:
%
%    parameters.sweep     sweep width in Hz
%
%    parameters.npoints   number of points for both dimensions
%
%    parameters.spins     nuclei on which the sequence runs,
%                         specified as {'1H'}, {'13C'}, etc.
%
%    parameters.tmix      mixing time, seconds
%
%    parameters.gamp      gradient amplitude, T/m
%
%    parameters.sal_ang   flip angle of the saltire chirp (degrees)
%
%    parameters.sal_dur   pulse width of saltire chirp (s)
%
%    parameters.sal_del   chirp pulse gradient duration (s)
%
%    parameters.sal_swp   sweep width of saltire chirp (Hz)
%
%    parameters.sal_npt   number of points in the saltire chirp
%
%    parameters.sal_smf   saltire chirp smoothing factor
%
%    H                    Fokker-Planck Hamiltonian, received
%                         from the imaging context
%
%    R                    Fokker-Planck relaxation superoperator,
%                         received from the imaging context
%
%    K                    Fokker-Planck kinetics superoperator,
%                         received from the imaging context
%
%    G                    Fokker-Planck gradient superoperators,
%                         received from the imaging context
%
%    F                    Fokker-Planck diffusion and flow super-
%                         operator, received from the context
%
% Outputs:
%
%    fid  - two-dimensional free induction decay
%
% a.m.kenwright@durham.ac.uk
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=psycosy.m>

function fid=psycosy(spin_system,parameters,H,R,K,G,F)

% Check consistency
grumble(spin_system,parameters,H,R,K,F,G)

% Coherent evolution timestep
parameters.delta=1/(4*parameters.sweep);

% Compose Liouvillian
L=H+F+1i*R+1i*K;

% Get pulse operators
Lp=operator(spin_system,'L+',parameters.spins{1});
Lx=kron(speye(prod(parameters.npts)),(Lp+Lp')/2);
Ly=kron(speye(prod(parameters.npts)),(Lp-Lp')/2i);

% Get chirp pulse qaveform
[Cx,Cy]=chirp_pulse(parameters.sal_npt,parameters.sal_dur,...
                       parameters.sal_swp,parameters.sal_smf,'saltire');

% Normalize chirp amplitude                   
norm_factor=max(Cx); Cx=Cx/norm_factor; Cy=Cy/norm_factor;
    
% Q factor for saltire flip angle  
q_beta=-(2*log(cosd(parameters.sal_ang)/2+1/2))/pi; 
        
% RF field strength of saltire chirp based on Q factor (Hz)
rfbeta=sqrt(parameters.sal_swp*q_beta/(2*pi*parameters.sal_dur));    

% Calibrate chirps
Cx=2*pi*rfbeta*Cx; Cy=2*pi*rfbeta*Cy;

% Apply the first pulse
rho_stack=step(spin_system,Lx,parameters.rho0,pi/2);

% Run the first half of the t1 evolution
rho_stack=evolution(spin_system,L,[],rho_stack,0.5/parameters.sweep,...
                    parameters.npoints(1)-1,'trajectory');

% Select "+1" coherence
rho_stack=coherence(spin_system,rho_stack,{{parameters.spins{1},+1}});

% Apply the hard 180 pulse
rho_stack=step(spin_system,Lx,rho_stack,pi);

% Select "-1" coherence
rho_stack=coherence(spin_system,rho_stack,{{parameters.spins{1},-1}});

% Apply the 1st pulse of the PSYCHE element
durations=ones(size(Cx))*parameters.sal_dur/numel(Cx);
rho_stack=shaped_pulse_xy(spin_system,L+parameters.gamp*G{1},...
                          {Lx,Ly},{Cx,+Cy},durations,rho_stack,'expv-pwc');

% Run the diffusion time evolution
rho_stack=evolution(spin_system,L+parameters.gamp*G{1},[],rho_stack,...
                    parameters.sal_del-2*parameters.sal_dur,1,'final');

% Select "0" coherence
rho_stack=coherence(spin_system,rho_stack,{{parameters.spins{1},0}});

% Apply the 2nd pulse of the PSYCHE element
rho_stack=shaped_pulse_xy(spin_system,L+parameters.gamp*G{1},...
                          {Lx,Ly},{Cx,-Cy},durations,rho_stack,'expv-pwc');
    
% Run the second half of the t1 evolution
rho_stack=evolution(spin_system,L,[],rho_stack,0.5/parameters.sweep,...
                    parameters.npoints(1)-1,'refocus');

% Select "+1" coherence
rho_stack=coherence(spin_system,rho_stack,{{parameters.spins{1},+1}});

% Run the mixing time under full Liouvillian
rho_stack=evolution(spin_system,L,[],rho_stack,parameters.tmix,1,'final');

% Select "+1" coherence
rho_stack=coherence(spin_system,rho_stack,{{parameters.spins{1},+1}});

% Apply the final pulse
rho_stack=step(spin_system,Lx,rho_stack,pi/2);
 
% Run the F2 evolution under full Liouvillian
fid=evolution(spin_system,L,parameters.coil,rho_stack,1/parameters.sweep,...
              parameters.npoints(2)-1,'observable');

end

% Consistency enforcement
function grumble(spin_system,parameters,H,R,K,F,G)
if ~ismember(spin_system.bas.formalism,{'sphten-liouv'})
    error('this function is only available for sphten-liouv formalism.');
end
if (~isnumeric(H))||(~isnumeric(R))||(~isnumeric(K))||(~isnumeric(F))||...
   (~ismatrix(H))||(~ismatrix(R))||(~ismatrix(K))||(~ismatrix(F))
    error('H, R, K and F must be matrices.');
end
if (~all(size(H)==size(R)))||(~all(size(R)==size(K)))||(~all(size(K)==size(F)))
    error('H, R, K and F matrices must have the same dimension.');
end
if ~iscell(G)
    error('the G must be a 1x3 cell array.');
end
if numel(G)~=3
    error('G must contain three gradient operators.');
end
if (~isnumeric(G{1}))||(~isnumeric(G{2}))||(~isnumeric(G{3}))
    error('gradient operators must be numeric matrices.');
end
if ~isfield(parameters,'gamp')
    error('the gradient amplitude should be specified in parameters.gamp variable.');
elseif (~isnumeric(parameters.gamp))||(~isreal(parameters.gamp))||...
       (numel(parameters.gamp)~=1)||(~isfinite(parameters.gamp))
    error('parameters.gamp must be a finite real scalar.');
end
if ~isfield(parameters,'spins')
    error('working spins should be specified in parameters.spins variable.');
elseif (~iscell(parameters.spins))||(numel(parameters.spins)~=1)||...
       (~ischar(parameters.spins{1}))
    error('parameters.spins must be a one-element cell array of character strings.');
end
if ~isfield(parameters,'rho0')
    error('initial state should be specified in parameters.rho0 variable.');
elseif (~isnumeric(parameters.rho0))||(~iscolumn(parameters.rho0))||...
       (size(parameters.rho0,1)~=size(H,1))
    error('parameters.rho0 must be a state vector of the same dimension as H.');
end
if ~isfield(parameters,'coil')
    error('detection state should be specified in parameters.coil variable.');
elseif (~isnumeric(parameters.coil))||(~iscolumn(parameters.coil))||...
       (size(parameters.coil,1)~=size(H,1))
    error('parameters.coil must be a state vector of the same dimension as H.');
end
if ~isfield(parameters,'sweep')
    error('sweep width should be specified in parameters.sweep variable.');
elseif (~isnumeric(parameters.sweep))||(~isreal(parameters.sweep))||...
       (numel(parameters.sweep)~=1)||(~isfinite(parameters.sweep))||(parameters.sweep<=0)
    error('parameters.sweep must be a positive real scalar.');
end
if ~isfield(parameters,'npoints')
    error('point counts should be specified in parameters.npoints variable.');
elseif (~isnumeric(parameters.npoints))||(~isreal(parameters.npoints))||...
       (numel(parameters.npoints)~=2)||any(~isfinite(parameters.npoints))||...
       any(parameters.npoints<2)||any(mod(parameters.npoints,1)~=0)
    error('parameters.npoints must contain two integer elements greater than 1.');
end
if ~isfield(parameters,'tmix')
    error('mixing time should be specified in parameters.tmix variable.');
elseif (~isnumeric(parameters.tmix))||(~isreal(parameters.tmix))||...
       (numel(parameters.tmix)~=1)||(~isfinite(parameters.tmix))||(parameters.tmix<0)
    error('parameters.tmix must be a non-negative real scalar.');
end
if ~isfield(parameters,'sal_ang')
    error('saltire flip angle should be specified in parameters.sal_ang variable.');
elseif (~isnumeric(parameters.sal_ang))||(~isreal(parameters.sal_ang))||...
       (numel(parameters.sal_ang)~=1)||(~isfinite(parameters.sal_ang))
    error('parameters.sal_ang must be a finite real scalar.');
end
if ~isfield(parameters,'sal_dur')
    error('saltire duration should be specified in parameters.sal_dur variable.');
elseif (~isnumeric(parameters.sal_dur))||(~isreal(parameters.sal_dur))||...
       (numel(parameters.sal_dur)~=1)||(~isfinite(parameters.sal_dur))||(parameters.sal_dur<=0)
    error('parameters.sal_dur must be a positive real scalar.');
end
if ~isfield(parameters,'sal_del')
    error('gradient duration should be specified in parameters.sal_del variable.');
elseif (~isnumeric(parameters.sal_del))||(~isreal(parameters.sal_del))||...
       (numel(parameters.sal_del)~=1)||(~isfinite(parameters.sal_del))||...
       (parameters.sal_del<=2*parameters.sal_dur)
    error('parameters.sal_del must be a real scalar greater than twice parameters.sal_dur.');
end
if ~isfield(parameters,'sal_swp')
    error('saltire sweep width should be specified in parameters.sal_swp variable.');
elseif (~isnumeric(parameters.sal_swp))||(~isreal(parameters.sal_swp))||...
       (numel(parameters.sal_swp)~=1)||(~isfinite(parameters.sal_swp))||(parameters.sal_swp<=0)
    error('parameters.sal_swp must be a positive real scalar.');
end
if ~isfield(parameters,'sal_npt')
    error('saltire point count should be specified in parameters.sal_npt variable.');
elseif (~isnumeric(parameters.sal_npt))||(~isreal(parameters.sal_npt))||...
       (numel(parameters.sal_npt)~=1)||(~isfinite(parameters.sal_npt))||...
       (parameters.sal_npt<1)||(mod(parameters.sal_npt,1)~=0)
    error('parameters.sal_npt must be a positive integer scalar.');
end
if ~isfield(parameters,'sal_smf')
    error('saltire smoothing factor should be specified in parameters.sal_smf variable.');
elseif (~isnumeric(parameters.sal_smf))||(~isreal(parameters.sal_smf))||...
       (numel(parameters.sal_smf)~=1)||(~isfinite(parameters.sal_smf))||...
       (parameters.sal_smf<0)||(parameters.sal_smf>50)
    error('parameters.sal_smf must be a real scalar between 0 and 50.');
end
end

% As a test of whether machines can think, the original 
% Turing test has been criticised on many counts, not the
% least being the ethical issue of how to treat human 
% beings who consistently fail to convince other human 
% beings of their humanity.
% 
% Seth Lloyd

