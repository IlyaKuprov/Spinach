% Two-pulse ESEEM echo; ideal pi/2 preparation and midpoint pi pulse.
% Syntax:
%
%                 E=cce_eseem(spin_system,parameters,H,R,K)
%
% parameters.times contains total sequence durations in seconds. transition
% selects the observed two-state transition; rho0/coil or preparation may
% override the prepared density and detector. normalized defaults to true
% and divides by the actual detected initial reference. nloops is required
% for CPMG. bath_pumped_hahn requires pump.center, halfwidth, units (Hz or
% rad/s), fraction and allow_empty; optional sites/species restrict the bath.
% Ideal refocusing acts on the selected transition about x. H,R,K are the
% conventional generators or explicitly tagged projected branch generators.
%
% ilya.kuprov@weizmann.ac.il

function E=cce_eseem(spin_system,parameters,H,R,K)

% Preparation and detection are owned by this experiment
grumble(spin_system,parameters,H,R,K);
if ~isfield(parameters,'normalized'), parameters.normalized=true; end
[rho0,coil,reference,meta]=cce_observation(spin_system,parameters);
if ~iscell(coil), coil={coil}; end
projected=isstruct(H); E=complex(zeros(numel(parameters.times),numel(reference)));
if projected, generator=H; else, generator=H+1i*R+1i*K; end

% Build ideal observed-spin transition pulses without a global bath basis
model=spin_system.model; dc=spin_system.observed_dimension;
db=spin_system.bath_dimension;
if ~isfield(model,'transition')
    error('cce_eseem:transition','An explicit transition is required for ideal refocusing.');
end
a=model.transition.alpha; b=model.transition.beta;
rotation=eye(dc)-a*a'-b*b'-1i*(a*b'+b*a');
rotation=kron(rotation,eye(db));
fractions=0.5;

% At zero elapsed time no control is applied, matching the reference signal
for n=1:numel(parameters.times)
    rho=rho0; previous=0; swapped=false;
    if parameters.times(n)>0
        for j=1:numel(fractions)
            segment=generator;
            if projected&&swapped, segment.left=H.right; segment.right=H.left; end
            duration=(fractions(j)-previous)*parameters.times(n);
            rho=cce_evolution(spin_system,segment,[],rho,duration,1,'final');
            if projected
                swapped=~swapped;
            else
                rho=reshape(rho,dc*db,dc*db); rho=rotation*rho*rotation';
                if strcmp(spin_system.formalism,'zeeman-liouv'), rho=rho(:); end
            end
            previous=fractions(j);
        end
    end
    segment=generator;
    if projected&&swapped, segment.left=H.right; segment.right=H.left; end
    rho=cce_evolution(spin_system,segment,[],rho,(1-previous)*parameters.times(n),1,'final');
    for k=1:numel(coil)
        if isvector(rho), value=coil{k}'*rho;
        else, value=trace(coil{k}'*rho); end
        if meta.normalized, value=value/reference(k); end
        E(n,k)=value;
    end
end

end

% Validate the registered sequence contract before preparation
function grumble(spin_system,parameters,H,R,K)
if ~isstruct(spin_system)||~isfield(spin_system,'kind')||...
   ~strcmp(spin_system.kind,'cce_subsystem')||~isstruct(parameters)||...
   ~isfield(parameters,'times')||~isnumeric(parameters.times)||...
   ~isreal(parameters.times)||~isvector(parameters.times)||...
   any(~isfinite(parameters.times))||any(parameters.times<0)||any(diff(parameters.times)<0)
    error('cce_eseem:input','A native retained subsystem and increasing nonnegative times are required.');
end
if isfield(parameters,'sequence')&&~isempty(parameters.sequence)
    error('cce_eseem:sequence','This named experiment defines its own controls.');
end
if isfield(parameters,'output')&&strcmp(parameters.output,'density')
    error('cce_eseem:density','This experiment returns detected signal or observables.');
end
if ~(isnumeric(H)||isstruct(H))||~isnumeric(R)||~isnumeric(K)
    error('cce_eseem:generators','Supply conventional or tagged projected generators.');
end

end


