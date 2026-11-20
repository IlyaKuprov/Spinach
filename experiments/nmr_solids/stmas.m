% Shifted-echo satellite-transition MAS with three finite x pulses.
% Call with the singlerot context and a single 87Rb spin in a
% complete zeeman-liouv basis. Syntax:
%
%             fid=stmas(spin_system,parameters,H,R,K)
%
% Parameters:
%
%     spin_system - Spinach system
%     parameters  - singlerot parameters, with rho0 as Lz and coil
%                   as CT-, pulse_amp (three amplitudes in rad/s),
%                   pulse_dur (three durations in seconds), tau
%                   (echo delay in seconds), and npoints=[ni np]
%     H,R,K       - rotor-space Hamiltonian, relaxation, and kinetics
%
% Outputs:
%
%     fid         - complex np-by-ni array, direct time in rows
%
% The first satellite evolution increment is one rotor period
% minus half of each of the first two pulse durations. The CT+
% delay is tau plus 8/9 of that free interval. Direct acquisition
% starts immediately after inversion and samples once per turn.
% The indirect sampling interval is 17/9 of a rotor period.
%
% talos@spindynamics.org

function fid=stmas(spin_system,parameters,H,R,K)

% Check consistency
grumble(spin_system,parameters,H,R,K);

% Set the centre-to-centre rotor clock and the two free intervals
period=1/abs(parameters.rate);
duration=parameters.pulse_dur;
sat_delay=period-(duration(1)+duration(2))/2;
ct_delay=parameters.tau+(8/9)*sat_delay;

% Compose the drift and lift spin operators and transition projectors
L=H+1i*R+1i*K;
space_unit=speye(parameters.spc_dim);
Lx=kron(space_unit,operator(spin_system,'Lx',1));
satellite=kron(space_unit,stmas_filter(spin_system,1,[1 2; 3 4]));
central=kron(space_unit,stmas_filter(spin_system,1,[2 3]));

% Apply excitation once and retain the positive-order satellites
rho_sat=step(spin_system,L+parameters.pulse_amp(1)*Lx,parameters.rho0,duration(1));
rho_sat=satellite*rho_sat;

% Propagate the first satellite interval and reuse one turn for later increments
rho_sat=step(spin_system,L,rho_sat,sat_delay);
sat_stack=evolution(spin_system,L,[],rho_sat,period,...
                    parameters.npoints(1)-1,'trajectory');
ct_stack=central*step(spin_system,L+parameters.pulse_amp(2)*Lx,sat_stack,duration(2));

% Split the common CT delay into bounded rotor-period steps
turns=floor(ct_delay/period);
ct_stack=evolution(spin_system,L,[],ct_stack,period,turns,'final');
ct_stack=step(spin_system,L,ct_stack,ct_delay-turns*period);

% Apply the incremental CT delay by powers of a short propagator
increment=propagator(spin_system,L,(8/9)*period);
indices=0:parameters.npoints(1)-1;
while any(indices)

    % Apply the current binary power to the relevant increments
    selected=logical(mod(indices,2));
    ct_stack(:,selected)=increment*ct_stack(:,selected);
    indices=floor(indices/2);
    if any(indices)
        increment=increment*increment;
    end

end

% Invert at the actual rotor phases and acquire every increment
rho=step(spin_system,L+parameters.pulse_amp(3)*Lx,ct_stack,duration(3));
fid=evolution(spin_system,L,parameters.coil,rho,period,...
              parameters.npoints(2)-1,'observable');

end

% Consistency enforcement
function grumble(spin_system,parameters,H,R,K)
if ~strcmp(spin_system.bas.formalism,'zeeman-liouv')||...
   ~isequal(spin_system.comp.isotopes,{'87Rb'})
    error('STMAS requires a zeeman-liouv basis for one 87Rb spin.');
end
if ~all(isfield(parameters,{'spc_dim','rho0','coil','rate','pulse_amp','pulse_dur','npoints','tau'}))
    error('STMAS context, timing, or pulse parameters are missing.');
end
if ~isnumeric(H)||~isnumeric(R)||~isnumeric(K)||...
   ~isequal(size(H),size(R),size(K),[16 16]*parameters.spc_dim)||...
   ~isequal(size(parameters.rho0),size(parameters.coil),[size(H,1) 1])
    error('drift and state dimensions must agree with the lifted spin basis.');
end
if ~isnumeric(parameters.pulse_amp)||~isreal(parameters.pulse_amp)||...
   ~isequal(size(parameters.pulse_amp),[1 3])||any(~isfinite(parameters.pulse_amp))
    error('pulse_amp must contain three finite real amplitudes.');
end
if ~isnumeric(parameters.pulse_dur)||~isreal(parameters.pulse_dur)||...
   ~isequal(size(parameters.pulse_dur),[1 3])||any(~isfinite(parameters.pulse_dur))||...
   any(parameters.pulse_dur<=0)
    error('pulse_dur must contain three finite positive durations.');
end
if ~isnumeric(parameters.npoints)||~isreal(parameters.npoints)||...
   ~isequal(size(parameters.npoints),[1 2])||...
   any(mod(parameters.npoints,1)~=0)||any(parameters.npoints<1)
    error('npoints must contain positive indirect and direct point counts.');
end
if ~isnumeric(parameters.rate)||~isreal(parameters.rate)||...
   ~isscalar(parameters.rate)||~isfinite(parameters.rate)||parameters.rate==0
    error('rate must be a finite nonzero spinning frequency.');
end
if ~isnumeric(parameters.tau)||~isreal(parameters.tau)||...
   ~isscalar(parameters.tau)||~isfinite(parameters.tau)||parameters.tau<0
    error('tau must be a finite nonnegative delay.');
end
if sum(parameters.pulse_dur(1:2))/2>1/abs(parameters.rate)
    error('pulse centres imply a negative first satellite free delay.');
end
end


