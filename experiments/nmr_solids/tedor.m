% Finite-pulse TEDOR for the Liouville-space singlerot context.
% The direct variant prepares the first channel and acquires the
% second with a stored rotor-period propagator. The two ZF variants
% prepare and detect the first channel; they differ in selection.
%
% Syntax:
%
%                    signal=tedor(spin_system,parameters,H,R,K)
%
% Parameters:
%
%    parameters.spins     - cell row of two distinct spin-1/2 isotope
%                           labels, prepared channel first
%    parameters.sequence  - 'direct', 'zf-cycle', or 'zf-filter'
%    parameters.rate      - non-zero signed MAS rate in Hz
%    parameters.ncycles   - non-negative integer row of mixing counts
%    parameters.rf        - positive common pulse amplitude in Hz
%    parameters.npoints   - positive integer direct acquisition length
%    parameters.rho0      - lifted starting column state
%    parameters.coil      - lifted detection column state
%    parameters.spc_dim   - spatial dimension supplied by singlerot
%    H,R,K                - explicit matrices from singlerot
%
% Outputs:
%
%    signal               - raw complex matrix; one column per count.
%                           Direct rows are rotor-period FID samples,
%                           starting immediately after conversion.
%                           ZF variants return one echo row.
%
% talos@spindynamics.org

function signal=tedor(spin_system,parameters,H,R,K)

% Check consistency
grumble(spin_system,parameters,H,R,K);

% Lift transverse controls and retain the drift during finite pulses
L=H+1i*R+1i*K;
space_id=speye(parameters.spc_dim);
Ix=kron(space_id,operator(spin_system,'Lx',parameters.spins{1}));
Sx=kron(space_id,operator(spin_system,'Lx',parameters.spins{2}));
Sy=kron(space_id,operator(spin_system,'Ly',parameters.spins{2}));
rf=2*pi*parameters.rf;
width=0.5/parameters.rf;
period=1/abs(parameters.rate);
half_gap=period/2-width;
signal=complex(zeros(1,numel(parameters.ncycles)));

% Direct transfer and rotor-synchronised acquisition
if strcmp(parameters.sequence,'direct')

    % Use separate receiver and sender inversion trains
    times=[half_gap/2 width half_gap width half_gap/2];
    mix={L,L+rf*Sx,L,L+rf*Sx,L};
    detect={L,L+rf*Ix,L,L+rf*Ix,L};
    signal=complex(zeros(parameters.npoints,numel(parameters.ncycles)));
    for n=1:numel(parameters.ncycles)

        % Prepare antiphase coherence, then convert on both channels
        rho=parameters.rho0;
        for k=1:parameters.ncycles(n)
            rho=prop_train(spin_system,mix,times,rho);
        end
        rho=step(spin_system,L+rf*(Ix+Sx),rho,width/2);

        % Acquire before applying each stored receiver-period block
        for k=1:parameters.npoints
            signal(k,n)=parameters.coil'*rho;
            if k<parameters.npoints
                rho=prop_train(spin_system,detect,times,rho);
            end
        end

    end
    return;

end

% Set up the two ZF recoupling arms and their central echo
first={L,L+rf*Sx,L,L+rf*Sy};
second={L+rf*Sx,L,L+rf*Sy,L};
first_times=[half_gap width half_gap width];
second_times=[width half_gap width half_gap];
echo={L,L+rf*Ix,L};
echo_times=[period-width/2 width period-width/2];

% Opposite and equal phase combinations of the two ideal conversion pulses
for n=1:numel(parameters.ncycles)

    % Prepare before conversion using the source order of both arms
    count=parameters.ncycles(n);
    rho=parameters.rho0;
    for k=1:count
        rho=prop_train(spin_system,first,first_times,rho);
    end
    rho=prop_train(spin_system,echo,echo_times,rho);
    for k=1:count
        rho=prop_train(spin_system,second,second_times,rho);
    end

    % The two ideal quarter turns commute and combine into pi rotations
    minus=step(spin_system,Ix,rho,pi);
    if strcmp(parameters.sequence,'zf-filter')

        % Select observed single-quantum and zero partner coherence
        rho=coherence(spin_system,minus,{{parameters.spins{1},[-1 1]},...
                                         {parameters.spins{2},0}});
    else

        % Add the equal-phase branch with positive receiver and subtract the other
        rho=step(spin_system,Ix+Sx,rho,pi)-minus;

    end

    % Repeat the same recoupling sandwich and detect the signed echo
    for k=1:count
        rho=prop_train(spin_system,first,first_times,rho);
    end
    rho=prop_train(spin_system,echo,echo_times,rho);
    for k=1:count
        rho=prop_train(spin_system,second,second_times,rho);
    end
    signal(n)=parameters.coil'*rho;

end

end

% Apply a chronological train of finite-duration generators
function rho=prop_train(spin_system,generators,times,rho)

% Advance the state through each pulse or free interval
for n=1:numel(times)
    rho=step(spin_system,generators{n},rho,times(n));
end

end

% Consistency enforcement
function grumble(spin_system,parameters,H,R,K)
if ~ismember(spin_system.bas.formalism,{'sphten-liouv','zeeman-liouv'})
    error('TEDOR requires sphten-liouv or zeeman-liouv formalism.');
end
if ~isnumeric(parameters.rho0)||~iscolumn(parameters.rho0)||...
   ~isnumeric(parameters.coil)||~isequal(size(parameters.coil),size(parameters.rho0))
    error('parameters.rho0 and parameters.coil must be matching column vectors.');
end
if ~isnumeric(H)||~isnumeric(R)||~isnumeric(K)||...
   ~isequal(size(H),size(R),size(K),[numel(parameters.rho0) numel(parameters.rho0)])
    error('H, R, and K must be explicit matrices matching the lifted state.');
end
if ~iscell(parameters.spins)||~isequal(size(parameters.spins),[1 2])||...
   ~all(cellfun(@ischar,parameters.spins))||...
   any(~ismember(parameters.spins,spin_system.comp.isotopes))||...
   strcmp(parameters.spins{1},parameters.spins{2})||...
   spin_system.comp.nspins~=2||any(spin_system.comp.mults~=2)
    error('TEDOR requires the two distinct spin-1/2 isotopes of the system.');
end
if ~ischar(parameters.sequence)||~ismember(parameters.sequence,{'direct','zf-cycle','zf-filter'})
    error('Unknown TEDOR sequence.');
end
validateattributes(parameters.rate,{'numeric'},{'scalar','real','finite','nonzero'});
validateattributes(parameters.rf,{'numeric'},{'scalar','real','finite','positive'});
validateattributes(parameters.ncycles,{'numeric'},{'row','nonempty','real','finite','nonnegative','integer'});
validateattributes(parameters.npoints,{'numeric'},{'scalar','real','finite','positive','integer'});
if parameters.rf<abs(parameters.rate)
    error('The pi pulse must fit within half a rotor period.');
end
end


