% Dipolar-shift (DIPSHIFT) echo for a heteronuclear pair under MAS,
% with one recoupling block in each arm. Use with the singlerot
% context in Liouville space. Syntax:
%
%             signal=dipshift(spin_system,parameters,H,R,K)
%
% Parameters:
%
%    parameters.spins       - {recoupled, observed} isotope labels
%
%    parameters.rate        - signed rotor rate, Hz
%
%    parameters.time_shift  - nonempty row of pulse displacements,
%                             in seconds, in [0,half a rotor period]
%
%    parameters.xy_block    - logical scalar; true for X,Y blocks
%                             in both arms, false for a single X
%                             pulse in each arm
%
%    parameters.pulse_width - pi pulse width in seconds, in
%                             [0,half a rotor period] for X,Y
%                             blocks; must be zero for X-only
%
%    parameters.rho0        - initial state, usually observed Lx
%
%    parameters.coil        - detection state, usually observed L+
%
%    H,R,K                  - Hamiltonian, relaxation, and kinetics
%                             generators, supplied by singlerot
%
% Outputs:
%
%    signal - complex echo row, one entry per time_shift value
%
% The echo drives X on the recoupled channel and Y on the observed
% channel simultaneously. A zero-width pulse retains its pi flip
% angle and acts at the same centre as the corresponding finite
% pulse. X-only and X,Y schedules have different pulse centres.
%
% talos@spindynamics.org

function signal=dipshift(spin_system,parameters,H,R,K)

% Check consistency
grumble(spin_system,parameters,H,R,K);

% Compose the drift and lift the RF generators into rotor space
L=H+1i*R+1i*K;
rotor_id=speye(parameters.spc_dim);
Ix=kron(rotor_id,operator(spin_system,'Lx',parameters.spins{1}));
Iy=kron(rotor_id,operator(spin_system,'Ly',parameters.spins{1}));
Sy=kron(rotor_id,operator(spin_system,'Ly',parameters.spins{2}));

% Put the pulse centres on their respective rotor schedules
period=1/abs(parameters.rate);
width=parameters.pulse_width;
if parameters.xy_block
    rf={Ix,Iy,Ix+Sy,Ix,Iy};
    centres=period*[1/2 1 3 5 11/2];
    shifts=[1 1 0 -1 -1];
    end_time=6*period;
else
    rf={Ix,Ix+Sy,Ix};
    centres=period*[1/2 2 7/2];
    shifts=[1 0 -1];
    end_time=4*period;
end

% Integrate pi rotations with drift, or use instantaneous rotations
if width>0

    % Scale the integrated generator to avoid dividing by a vanishing width
    for n=1:numel(rf)
        rf{n}=rf{n}+(width/pi)*L;
    end

end

% Restart each increment at the same rotor phase and prepared state
signal=complex(zeros(size(parameters.time_shift)));
for n=1:numel(signal)

    % Compute drift intervals between pulse edges
    times=centres+shifts*parameters.time_shift(n);
    delays=[times(1)-width/2 diff(times)-width end_time-times(end)-width/2];

    % Remove round-off at coincident pulse edges
    delays=max(delays,0);

    % Alternate drift and RF before detecting the echo
    rho=parameters.rho0;
    for k=1:numel(rf)
        rho=step(spin_system,L,rho,delays(k));
        rho=step(spin_system,rf{k},rho,pi);
    end
    rho=step(spin_system,L,rho,delays(end));
    signal(n)=parameters.coil'*rho;

end

end

% Consistency enforcement
function grumble(spin_system,parameters,H,R,K)
if ~ismember(spin_system.bas.formalism,{'sphten-liouv','zeeman-liouv'})
    error('dipshift requires a Liouville-space formalism.');
end
if ~isnumeric(H)||~isnumeric(R)||~isnumeric(K)||...
   ~ismatrix(H)||~ismatrix(R)||~ismatrix(K)||...
   ~isequal(size(H),size(R),size(K),[size(H,1) size(H,1)])
    error('H, R, and K must be square operators of the same size.');
end
if ~iscell(parameters.spins)||~isequal(size(parameters.spins),[1 2])||...
   ~all(cellfun(@ischar,parameters.spins))||...
   ~all(ismember(parameters.spins,spin_system.comp.isotopes))||...
   strcmp(parameters.spins{1},parameters.spins{2})
    error('parameters.spins must contain two distinct isotope labels present in the system.');
end
validateattributes(parameters.rate,{'numeric'},{'real','finite','scalar','nonzero'});
half_period=0.5/abs(parameters.rate);
validateattributes(parameters.time_shift,{'numeric'},...
                   {'real','finite','row','nonempty','>=',0,'<=',half_period});
validateattributes(parameters.xy_block,{'logical'},{'scalar'});
validateattributes(parameters.pulse_width,{'numeric'},...
                   {'real','finite','scalar','>=',0,'<=',half_period});
if ~parameters.xy_block&&(parameters.pulse_width~=0)
    error('the X-only schedule requires zero pulse width.');
end
validateattributes(parameters.rho0,{'numeric'},{'column','numel',size(H,1)});
validateattributes(parameters.coil,{'numeric'},{'column','numel',size(H,1)});
end


