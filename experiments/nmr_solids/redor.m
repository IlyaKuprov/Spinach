% Rotational-echo double-resonance (REDOR) experiment with ideal
% hard pi pulses. The observed channel is refocused once per rotor
% period, and the dephasing channel is pulsed at rotor-period
% boundaries. The sequence reports the full rotational echo, the
% dephased echo, and their difference as a function of the number of
% rotor cycles. To be called from singlerot context. Further
% information in:
%
%                 https://doi.org/10.1016/0022-2364(89)90280-1
%                 https://doi.org/10.1006/jmre.2000.2128
%
% Syntax:
%
%          redor_curve=redor(spin_system,parameters,H,R,K)
%
% Parameters:
%
%    parameters.spins        - observed and dephasing spins,
%                              e.g. {'13C','15N'} for 13C{15N}
%                              REDOR
%
%    parameters.ncycles      - row vector with numbers of rotor
%                              cycles in the REDOR evolution time
%
%    parameters.rate         - MAS rate in Hz
%
%    parameters.rho0         - initial state vector, usually
%                              transverse magnetisation on the
%                              observed spin
%
%    parameters.coil         - detection state vector, usually on
%                              the observed spin
%
%    parameters.decouple     - nuclei to decouple analytically
%                              during REDOR evolution, supplied as
%                              a cell array of isotope strings
%
%    parameters.refocus_phase - phase of the observed-channel pi
%                              refocusing pulses, radians; defaults
%                              to 0
%
%    parameters.pulse_phase  - phases of the dephasing-channel pi
%                              pulses, radians; the vector is cycled
%                              through the pulse train, and defaults
%                              to [0 pi/2]
%
%    H  - Hamiltonian matrix, received from context function
%
%    R  - relaxation superoperator, received from context function
%
%    K  - kinetics superoperator, received from context function
%
% Outputs:
%
%    redor_curve(1,:)        - full echo S0, with observed-channel
%                              refocusing pulses only
%
%    redor_curve(2,:)        - dephased echo S, with observed-channel
%                              refocusing and dephasing-channel pi
%                              pulses
%
%    redor_curve(3,:)        - REDOR difference, S0-S
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=redor.m>

function redor_curve=redor(spin_system,parameters,H,R,K)

% Set default decoupling
if ~isfield(parameters,'decouple')
    parameters.decouple={};
end

% Set default refocusing phase
if ~isfield(parameters,'refocus_phase')
    parameters.refocus_phase=0;
end

% Set default dephasing pulse phases
if ~isfield(parameters,'pulse_phase')
    parameters.pulse_phase=[0 pi/2];
end

% Check consistency
grumble(spin_system,parameters,H,R,K);

% Compose the Liouvillian
L=H+1i*R+1i*K;

% Apply analytical decoupling if requested
[L,rho0]=decouple(spin_system,L,parameters.rho0,parameters.decouple);

% Get observed-channel refocusing operator
Ip=operator(spin_system,'L+',parameters.spins{1});
Ix=(Ip+Ip')/2; Iy=(Ip-Ip')/2i;
ref_oper=cos(parameters.refocus_phase)*Ix+...
         sin(parameters.refocus_phase)*Iy;
ref_oper=kron(speye(parameters.spc_dim),ref_oper);

% Get dephasing-channel pulse operators
Ip=operator(spin_system,'L+',parameters.spins{2});
Ix=(Ip+Ip')/2; Iy=(Ip-Ip')/2i;
Ix=kron(speye(parameters.spc_dim),Ix);
Iy=kron(speye(parameters.spc_dim),Iy);

% Precompute phased dephasing-channel pulse operators
deph_opers=cell(size(parameters.pulse_phase));
for n=1:numel(parameters.pulse_phase)
    deph_opers{n}=cos(parameters.pulse_phase(n))*Ix+...
                  sin(parameters.pulse_phase(n))*Iy;
end

% Get rotor-synchronised timing
rotor_period=1/abs(parameters.rate);
half_period=rotor_period/2;

% Preallocate echo arrays
s0_echo=zeros(1,max(parameters.ncycles)+1);
s_echo=zeros(1,max(parameters.ncycles)+1);

% Store the zero-time point
s0_echo(1)=parameters.coil'*rho0;
s_echo(1)=parameters.coil'*rho0;

% Initialise reference and dephased trajectories
rho_s0=rho0; rho_s=rho0; pulse_count=0;

% Step through the requested REDOR evolution window
for n=1:max(parameters.ncycles)

    % Propagate both echoes over the first half-period
    rho_s0=step(spin_system,L,rho_s0,half_period);
    rho_s=step(spin_system,L,rho_s,half_period);

    % Refocus the observed channel
    rho_s0=step(spin_system,ref_oper,rho_s0,pi);
    rho_s=step(spin_system,ref_oper,rho_s,pi);

    % Propagate both echoes over the second half-period
    rho_s0=step(spin_system,L,rho_s0,half_period);
    rho_s=step(spin_system,L,rho_s,half_period);

    % Apply the dephasing-channel pi pulse to the S echo
    pulse_count=pulse_count+1;
    pulse_idx=mod(pulse_count-1,numel(deph_opers))+1;
    rho_s=step(spin_system,deph_opers{pulse_idx},rho_s,pi);

    % Store the rotor-synchronised echo intensities
    s0_echo(n+1)=parameters.coil'*rho_s0;
    s_echo(n+1)=parameters.coil'*rho_s;

end

% Pick out the requested rotor-cycle counts
s0_echo=s0_echo(parameters.ncycles+1);
s_echo=s_echo(parameters.ncycles+1);

% Assemble the REDOR curve
redor_curve=[s0_echo; s_echo; s0_echo-s_echo];

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
if ~isfield(parameters,'spins')
    error('working spins should be specified in parameters.spins variable.');
end
if (~iscell(parameters.spins))||(numel(parameters.spins)~=2)
    error('parameters.spins cell array should have exactly two elements.');
end
if any(~cellfun(@ischar,parameters.spins))
    error('elements of parameters.spins cell array must be strings.');
end
if any(~ismember(parameters.spins,spin_system.comp.isotopes))
    error('parameters.spins contains isotopes that are not present in the system.');
end
if ~isfield(parameters,'ncycles')
    error('numbers of rotor cycles should be specified in parameters.ncycles variable.');
end
if (~isnumeric(parameters.ncycles))||(~isreal(parameters.ncycles))||...
   (~isrow(parameters.ncycles))||isempty(parameters.ncycles)||...
   any(parameters.ncycles<0)||any(mod(parameters.ncycles,1)~=0)
    error('parameters.ncycles must be a row vector of non-negative integers.');
end
if ~isfield(parameters,'rate')
    error('MAS rate should be specified in parameters.rate variable.');
end
if (~isnumeric(parameters.rate))||(~isreal(parameters.rate))||...
   (~isscalar(parameters.rate))||(parameters.rate==0)
    error('parameters.rate must be a non-zero real scalar.');
end
if ~isfield(parameters,'rho0')
    error('initial state must be specified in parameters.rho0 variable.');
end
if (~isnumeric(parameters.rho0))||(~iscolumn(parameters.rho0))
    error('parameters.rho0 must be a state vector.');
end
if ~isfield(parameters,'coil')
    error('detection state must be specified in parameters.coil variable.');
end
if (~isnumeric(parameters.coil))||(~iscolumn(parameters.coil))
    error('parameters.coil must be a state vector.');
end
if numel(parameters.rho0)~=numel(parameters.coil)
    error('parameters.rho0 and parameters.coil must have the same dimension.');
end
if ~isfield(parameters,'spc_dim')
    error('spatial problem dimension must be specified in parameters.spc_dim variable.');
end
if (~isnumeric(parameters.spc_dim))||(~isreal(parameters.spc_dim))||...
   (~isscalar(parameters.spc_dim))||(parameters.spc_dim<1)||...
   (mod(parameters.spc_dim,1)~=0)
    error('parameters.spc_dim must be a positive integer.');
end
if ~iscell(parameters.decouple)
    error('parameters.decouple must be a cell array of strings.');
end
if any(~cellfun(@ischar,parameters.decouple))
    error('elements of parameters.decouple cell array must be strings.');
end
if any(~ismember(parameters.decouple,spin_system.comp.isotopes))
    error('parameters.decouple contains isotopes that are not present in the system.');
end
if any(ismember(parameters.spins,parameters.decouple))
    error('parameters.decouple must not include REDOR working spins.');
end
if (~isnumeric(parameters.refocus_phase))||(~isreal(parameters.refocus_phase))||...
   (~isscalar(parameters.refocus_phase))
    error('parameters.refocus_phase must be a real scalar.');
end
if (~isnumeric(parameters.pulse_phase))||(~isreal(parameters.pulse_phase))||...
   (~isrow(parameters.pulse_phase))||isempty(parameters.pulse_phase)
    error('parameters.pulse_phase must be a non-empty row vector of real numbers.');
end
end

% The empires of the future are the empires of the mind.
%
% Winston Churchill

