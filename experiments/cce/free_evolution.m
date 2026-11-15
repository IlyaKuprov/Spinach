% Free evolution of a retained cluster with experiment-owned detection.
% Syntax:
%
%           answer=free_evolution(spin_system,parameters,H,R,K)
%
% parameters.times gives elapsed times in seconds. rho0 and coil are explicit
% observed-spin matrices, or transition selects preparation and detection.
% preparation may contain state and detector descriptors. normalized defaults
% to true for signals. output is signal, observable or density; density output
% is the column-vectorised observed density, one row per elapsed time.
% H,R,K are conventional angular-frequency generators or explicitly tagged
% projected branches from cce_generators. No pulse or receiver is implicit.
%
% ilya.kuprov@weizmann.ac.il

function answer=free_evolution(spin_system,parameters,H,R,K)

% Preparation and detection belong to this experiment, not to its context
grumble(spin_system,parameters,H,R,K);
if ~isfield(parameters,'normalized'), parameters.normalized=true; end
[rho0,coil,reference,meta]=cce_observation(spin_system,parameters);
if isstruct(H), generator=H; else, generator=H+1i*R+1i*K; end
density=isfield(parameters,'output')&&strcmp(parameters.output,'density');
if density
    if isstruct(H), error('free_evolution:density','Projected branches do not return a density matrix.'); end
    answer=complex(zeros(numel(parameters.times),spin_system.observed_dimension^2));
else
    answer=complex(zeros(numel(parameters.times),numel(reference)));
end

% Every requested elapsed time starts from the same prepared density matrix
for n=1:numel(parameters.times)
    rho=cce_evolution(spin_system,generator,[],rho0,parameters.times(n),1,'final');
    if density
        dc=spin_system.observed_dimension; db=spin_system.bath_dimension;
        rho=reshape(rho,dc*db,dc*db); observed=zeros(dc);
        for a=1:dc
            for b=1:dc
                observed(a,b)=trace(rho((a-1)*db+(1:db),(b-1)*db+(1:db)));
            end
        end
        answer(n,:)=observed(:).';
    else
        if ~iscell(coil), detectors={coil}; else, detectors=coil; end
        for k=1:numel(detectors)
            if isvector(rho), value=detectors{k}'*rho;
            else, value=trace(detectors{k}'*rho); end
            if meta.normalized, value=value/reference(k); end
            answer(n,k)=value;
        end
    end
end

end

% Validate the five-argument experiment contract and elapsed-time array
function grumble(spin_system,parameters,H,R,K)
if ~isstruct(spin_system)||~isfield(spin_system,'kind')||...
   ~strcmp(spin_system.kind,'cce_subsystem')||~isstruct(parameters)||...
   ~isfield(parameters,'times')||~isnumeric(parameters.times)||...
   ~isreal(parameters.times)||~isvector(parameters.times)||...
   any(~isfinite(parameters.times))||any(parameters.times<0)||any(diff(parameters.times)<0)
    error('free_evolution:input','A native retained subsystem and increasing nonnegative times are required.');
end
if ~(isnumeric(H)||isstruct(H))||~isnumeric(R)||~isnumeric(K)
    error('free_evolution:generators','Supply conventional or tagged projected generators.');
end
end


