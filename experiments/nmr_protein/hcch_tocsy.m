% HCCH-TOCSY pulse sequence from Figure 7.28 of Protein NMR Spectroscopy
% (2nd edition) using the bidirectional propagation method described in
%
%              http://dx.doi.org/10.1016/j.jmr.2014.04.002
%
% The sequence is hard-wired to work on 1H,13C proteins and uses PDB la-
% bels to select spins that will be affected by otherwise ideal pulses.
% F1 is 1H, F2 is C, F3 is H. Syntax:
%
%              fid=hcch_tocsy(spin_system,parameters,H,R,K)
%
% Parameters:
%
%    parameters.npoints     - a vector of three integers giving the
%                             number of points in the three temporal
%                             dimensions, ordered as [t1 t2 t3]
%
%    parameters.sweep       - a vector of three real numbers giving
%                             the sweep widths in the three frequen-
%                             cy dimensions, ordered as [f1 f2 f3]
%
%    parameters.J_ch        - 1H-13C J-coupling in Hz to be used for 
%                             magnetisation transfer, typically 140
%
%    parameters.delta       - coherence transfer evolution delay in
%                             seconds, see the pulse sequence diag-
%                             ram, typically 1.1e-3
%
%    parameters.lamp        - 1H-13C spin-lock amplitude in Hz, ty-
%                             pically 10000
%
%    parameters.sl_tmix     - spin-lock mixing time, typically 2e-3
%                             seconds
%
%    parameters.dipsi_dur   - DIPSI period duration, typically 
%                             22.5e-3 seconds
%
%    parameters.decouple_f3 - list of spins to be decoupled during
%                             the detection period, typically {'13C'}
%
%    H   - Hamiltonian matrix, received from context function
%
%    R   - relaxation superoperator, received from context function
%
%    K   - kinetics superoperator, received from context function
%
% Outputs:
%
%    fid - three-dimensional free induction decay
%
% Note: spin labels must be set to PDB atom IDs ('CA', 'HA', etc.) in
%       sys.labels for this sequence to work properly.
%
% TODO: whoever understands how phase cycles and quadratures work in
%       3D NMR is welcome to add a phase-sensitive version.
%
% m.walker@soton.ac.uk
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=hcch_tocsy.m>

function fid=hcch_tocsy(spin_system,parameters,H,R,K)

% Consistency check
grumble(spin_system,parameters,H,R,K);

% Compose Liouvillian
L=H+1i*R+1i*K;

% Coherent evolution timesteps
t1.nsteps=parameters.npoints(1); t1.timestep=1./parameters.sweep(1);
t2.nsteps=parameters.npoints(2); t2.timestep=1./parameters.sweep(2);
t3.nsteps=parameters.npoints(3); t3.timestep=1./parameters.sweep(3);

% J-coupling evolution time
tau_ch=abs(1/4*parameters.J_ch);

% Initial condition
rho0=state(spin_system,'Lz','1H','cheap');

% Detection state
coil=state(spin_system,'L+','1H','cheap');

% Pulse operators all protons
Hp=operator(spin_system,'L+','1H');
Hx=(Hp+Hp')/2; Hy=(Hp-Hp')/2i; 

% Pulse operators all carbons
Cp=operator(spin_system,'L+','13C');
Cx=(Cp+Cp')/2; Cy=(Cp-Cp')/2i;

% Selective pulse on CO carbons
COs=strcmp('C',spin_system.comp.labels); 
COp=operator(spin_system,'L+',find(COs));
COx=(COp+COp')/2;

%% Forward sim from rho0 up to t2 period 

% Pulse on 1H 
rho=step(spin_system,Hx,rho0,pi/2);

% Coherence selection on protons
rho=coherence(spin_system,rho,{{'1H',+1}});

% tau evolution
rho=evolution(spin_system,L,[],rho,tau_ch,1,'final');

% First part of t1 evolution
rho_stack=evolution(spin_system,L,[],rho,t1.timestep/2,t1.nsteps-1,'trajectory');

% Inversion pulse on 13C
rho_stack=step(spin_system,Cx,rho_stack,pi);

% The rest of t1 evolution
rho_stack=evolution(spin_system,L,[],rho_stack,t1.timestep/2,t1.nsteps-1,'refocus');    

% Inversion pulse on 1H
rho_stack=step(spin_system,Hx,rho_stack,pi);

% tau evolution
rho_stack=evolution(spin_system,L,[],rho_stack,tau_ch,1,'final');

% Pulses on 1H and 13C
rho_stack=step(spin_system,Hx+Cx,rho_stack,pi/2);

%% Backward sim from coil up to t2 period 

% Get decoupled evolution generator
[L_dec,~]=decouple(spin_system,L,[],parameters.decouple_f3);

% Backward detection on 1H
coil_stack=evolution(spin_system,L_dec,[],coil,-t3.timestep,...
                     t3.nsteps-1,'trajectory');
                 
% Select single quantum coherence
coil_stack=coherence(spin_system,coil_stack,{{'1H',1}});                 

% Backward tau evolution
coil_stack=evolution(spin_system,L,[],coil_stack,-tau_ch,1,'final');

% Backward inversion pulses on 1H and 13C
coil_stack=step(spin_system,Hx+Cx,coil_stack,-pi);

% Backward tau evolution
coil_stack=evolution(spin_system,L,[],coil_stack,-tau_ch,1,'final');

% Backward pulses on 1H and 13C
coil_stack=step(spin_system,Hx+Cx,coil_stack,-pi/2);

% Backward delta evolution
coil_stack=evolution(spin_system,L,[],coil_stack,-parameters.delta,1,'final');

% Backward inversion pulses on 1H and 13C
coil_stack=step(spin_system,Hx+Cx,coil_stack,-pi);

% Backward delta evolution
coil_stack=evolution(spin_system,L,[],coil_stack,-parameters.delta,1,'final');

% Effective isotropic mixing Liouvillian
spin_system=dictum(spin_system,{'1H'},'ignore');       % Ignore Zeeman term for 1H
spin_system=dictum(spin_system,{'13C'},'ignore');      % Ignore Zeeman term for 13C
spin_system=dictum(spin_system,{'1H','13C'},'strong'); % Enforce strong 1H-13C coupling
L_isomix=hamiltonian(spin_system)+1i*R+1i*K;       

% Backward evolution during isotropic mixing
coil_stack=evolution(spin_system,L_isomix,[],coil_stack,-parameters.dipsi_dur,1,'final');

% Backward spin lock
coil_stack=evolution(spin_system,L+2*pi*parameters.lamp*(Hy+Cy),[],coil_stack,...
                     -parameters.sl_tmix,1,'final');

% Backward delta evolution
coil_stack=evolution(spin_system,L,[],coil_stack,-parameters.delta,1,'final');

% Backward inversion pulse on 13C
coil_stack=step(spin_system,Cx,coil_stack,-pi);

%% Stitch the halves

% Coherence selection
rho_stack=coherence(spin_system,rho_stack,{{'13C',+1}});
coil_stack=coherence(spin_system,coil_stack,{{'13C',+1}});

% Bidirectional evolution and stitching
fid=stitch(spin_system,L,rho_stack,coil_stack,{COx,L,Hx},{pi,parameters.delta,pi},t1,t2,t3);

end

% Consistency enforcement
function grumble(spin_system,parameters,H,R,K)
if ~ismember(spin_system.bas.formalism,{'sphten-liouv'})
    error('this function is only available for sphten-liouv formalism.');
end
if (~isnumeric(H))||(~isnumeric(R))||(~isnumeric(K))||...
   (~ismatrix(H))||(~ismatrix(R))||(~ismatrix(K))
    error('H, R and K arguments must be matrices.');
end
if (~all(size(H)==size(R)))||(~all(size(R)==size(K)))
    error('H, R and K matrices must have the same dimension.');
end
if ~isfield(parameters,'sweep')
    error('sweep width should be specified in parameters.sweep variable.');
elseif numel(parameters.sweep)~=3
    error('parameters.sweep array should have exactly three elements.');
end
if ~isfield(parameters,'spins')
    error('working spins should be specified in parameters.spins variable.');
end
if (~isnumeric(parameters.sweep))||(~isvector(parameters.sweep))||...
   (~isreal(parameters.sweep))||(numel(parameters.sweep)~=3)
    error('parameters.sweep must be a vector of three real numbers.');
end
if ~isfield(parameters,'npoints')
    error('number of points should be specified in parameters.npoints variable.');
end
if (~isnumeric(parameters.npoints))||(~isvector(parameters.npoints))||...
   (~isreal(parameters.npoints))||(numel(parameters.npoints)~=3)||...
    any(parameters.npoints<1)||any(mod(parameters.npoints,1)~=0)
    error('parameters.npoints must be a vector of three positive integers.');
end
if ~isfield(parameters,'J_ch')
    error('scalar coupling should be specified in parameters.J_ch variable.');
elseif numel(parameters.J_ch)~=1
    error('parameters.J_ch array should have exactly one element.');
end
if ~isfield(parameters,'delta')
    error('delta delay should be specified in parameters.delta variable.');
elseif numel(parameters.delta)~=1
    error('parameters.delta array should have exactly one element.');
end
if ~isfield(parameters,'sl_tmix')
    error('spin lock duration should be specified in parameters.sl_tmix variable.');
elseif numel(parameters.sl_tmix)~=1
    error('parameters.sl_tmix array should have exactly one element.');
end
if ~isfield(parameters,'lamp')
    error('amplitude should be specified in parameters.lamp variable.');
elseif numel(parameters.lamp)~=1
    error('parameters.lamp array should have exactly one element.');
end
if ~isfield(parameters,'dipsi_dur')
    error('DIPSI duration should be specified in parameters.dipsi_dur variable.');
elseif numel(parameters.dipsi_dur)~=1
    error('parameters.dipsi_dur array should have exactly one element.');
end
end

% Talk is cheap, show me the code!
%
% Linus Torvalds

