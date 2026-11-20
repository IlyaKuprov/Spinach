% Repeated double-quantum-selected recoupling in the singlerot context.
% Compiles the finite RF schedule once, applies the requested number
% of blocks, selects total order +/-2, and repeats the same blocks
% before detection. Zero blocks return the initial overlap.
%
% Syntax: signal=relay_echo(spin_system,parameters,H,R,K)
%
% Parameters:
%    isotope       - irradiated isotope label
%    phases        - chronological pulse phases (radians)
%    widths        - two alternating widths (SR26), or one per phase
%    rf            - RF amplitude (rad/s)
%    sr26          - logical selecting the two-width compiler
%    counts        - row of successive counts starting at zero
%    spc_dim       - rotor-space dimension from singlerot
%    rho0,coil     - lifted preparation and receiver columns
%    H,R,K         - lifted Hamiltonian, relaxation, and kinetics
%
% Outputs: complex row in count order.
%
% talos@spindynamics.org

function signal=relay_echo(spin_system,parameters,H,R,K)

% Check consistency
grumble(spin_system,parameters,H,R,K);

% Form spatially replicated RF controls and the MAS drift
L=H+1i*R+1i*K;
space_id=speye(parameters.spc_dim);
Sx=kron(space_id,operator(spin_system,'Lx',parameters.isotope));
Sy=kron(space_id,operator(spin_system,'Ly',parameters.isotope));
if parameters.sr26

    % Reuse the distinct phase-duration segments of SR26
    [blocks,order]=rseq_compiler(spin_system,L,Sx,Sy,parameters.phases,...
                                 parameters.rf,parameters.widths,'90270_pulse');

else

    % Compile the S3 composite elements without changing their order
    blocks=cell(size(parameters.phases));
    order=1:numel(blocks);
    for n=order
        phase=parameters.phases(n);
        drive=parameters.rf*(cos(phase)*Sx+sin(phase)*Sy);
        blocks{n}=propagator(spin_system,L+drive,parameters.widths(n));
    end

end

% Accumulate the chronological product once
U=speye(size(L));
for n=1:numel(order)
    U=blocks{order(n)}*U;
end

% Advance preparation and receiver independently and select DQ coherence
rho=parameters.rho0;
coil=parameters.coil;
signal=zeros(size(parameters.counts));
signal(1)=coil'*rho;
for n=2:numel(signal)
    rho=U*rho;
    coil=U'*coil;
    signal(n)=coil'*coherence(spin_system,rho,{{'all',[-2 2]}});
end

end

% Consistency enforcement
function grumble(spin_system,parameters,H,R,K)
if ~ismember(spin_system.bas.formalism,{'sphten-liouv','zeeman-liouv'})
    error('DQ relay requires Liouville-space singlerot.');
end
if ~ischar(parameters.isotope)||~ismember(parameters.isotope,spin_system.comp.isotopes)
    error('Supply an irradiated isotope present in the system.');
end
if ~isnumeric(H)||~isnumeric(R)||~isnumeric(K)||...
   ~isequal(size(H),size(R),size(K),[numel(parameters.rho0) numel(parameters.rho0)])
    error('Lifted generators must match the state dimension.');
end
if ~iscolumn(parameters.rho0)||~isequal(size(parameters.coil),size(parameters.rho0))
    error('Preparation and receiver must be matching column states.');
end
validateattributes(parameters.phases,{'numeric'},{'column','real','finite','nonempty'});
validateattributes(parameters.widths,{'numeric'},{'row','real','finite','positive'});
validateattributes(parameters.rf,{'numeric'},{'scalar','real','finite','positive'});
validateattributes(parameters.sr26,{'logical'},{'scalar'});
if ~isequal(parameters.counts,0:numel(parameters.counts)-1)||isempty(parameters.counts)
    error('Supply successive counts from zero in a nonempty row.');
end
if parameters.sr26
    if numel(parameters.widths)~=2||mod(numel(parameters.phases),2)~=0
        error('The SR26 compiler needs paired phases and two alternating widths.');
    end
elseif numel(parameters.widths)~=numel(parameters.phases)
    error('Supply one duration for each phase.');
end
if parameters.sr26
    duration=sum(parameters.widths)*numel(parameters.phases)/2;
else
    duration=sum(parameters.widths);
end
rotations=duration*abs(parameters.rate);
if abs(rotations-round(rotations))>eps(rotations)*numel(parameters.phases)
    error('The repeated block duration must be an integer number of rotor periods.');
end
end


