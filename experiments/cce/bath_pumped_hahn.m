% Bath-pumped Hahn echo with an explicitly selected spectral window.
% Syntax:
%
%                 E=bath_pumped_hahn(spin_system,parameters,H,R,K)
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

function E=bath_pumped_hahn(spin_system,parameters,H,R,K)

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
    error('bath_pumped_hahn:transition','An explicit transition is required for ideal refocusing.');
end
a=model.transition.alpha; b=model.transition.beta;
rotation=eye(dc)-a*a'-b*b'-1i*(a*b'+b*a');
rotation=kron(rotation,eye(db));

% Select against canonical laboratory-frame bath frequencies in rad/s
spec=parameters.pump; ids=model.bath_ids;
if isfield(spec,'sites'), ids=intersect(ids,spec.sites,'stable'); end
if isfield(spec,'species'), ids=ids(ismember(model.isotopes(ids),cellstr(spec.species))); end
frequencies=zeros(size(ids));
for k=1:numel(ids)
    tensor=model.zeeman{ids(k)};
    if norm(tensor(3,1:2))>3*eps*norm(tensor(3,:))
        error('bath_pumped_hahn:axial','Spectral selection requires an axial canonical Zeeman field.');
    end
    frequencies(k)=tensor(3,3);
end
factor=1; if strcmp(spec.units,'Hz'), factor=2*pi; end
ids=ids(abs(frequencies-factor*spec.center)<=factor*spec.halfwidth);
if isempty(ids)&&~spec.allow_empty
    error('bath_pumped_hahn:empty','Pump window selected no bath spins.');
end
if isfield(parameters.cce,'meanfield')&&~strcmp(parameters.cce.meanfield,'off')&&...
   ~all(ismember(ids,spin_system.sites))
    error('bath_pumped_hahn:exterior','Pumping conditioned exterior spins requires the ensemble contract.');
end
ids=intersect(ids,spin_system.sites,'stable'); pump_operator=zeros(dc*db);
for id=reshape(ids,1,[])
    pump_operator=pump_operator+cce_operator(spin_system,'Lx',id);
end
if projected, pump_operator=pump_operator(1:db,1:db); end
pump_rotation=expm(-1i*pi*pump_operator);
[fractions,order]=sort([0.5 spec.fraction]); pulse_kind=[1 2]; pulse_kind=pulse_kind(order);

% At zero elapsed time no control is applied, matching the reference signal
for n=1:numel(parameters.times)
    rho=rho0; previous=0; swapped=false;
    if parameters.times(n)>0
        for j=1:numel(fractions)
            segment=generator;
            if projected&&swapped, segment.left=H.right; segment.right=H.left; end
            duration=(fractions(j)-previous)*parameters.times(n);
            rho=cce_evolution(spin_system,segment,[],rho,duration,1,'final');
            if pulse_kind(j)==1
                if projected
                    swapped=~swapped;
                else
                    rho=reshape(rho,dc*db,dc*db); rho=rotation*rho*rotation';
                    if strcmp(spin_system.formalism,'zeeman-liouv'), rho=rho(:); end
                end
            else
                if projected
                    rho=reshape(rho,db,db); rho=pump_rotation*rho*pump_rotation'; rho=rho(:);
                else
                    rho=reshape(rho,dc*db,dc*db); rho=pump_rotation*rho*pump_rotation';
                    if strcmp(spin_system.formalism,'zeeman-liouv'), rho=rho(:); end
                end
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
    error('bath_pumped_hahn:input','A native retained subsystem and increasing nonnegative times are required.');
end
if isfield(parameters,'sequence')&&~isempty(parameters.sequence)
    error('bath_pumped_hahn:sequence','This named experiment defines its own controls.');
end
if isfield(parameters,'output')&&strcmp(parameters.output,'density')
    error('bath_pumped_hahn:density','This experiment returns detected signal or observables.');
end
if ~(isnumeric(H)||isstruct(H))||~isnumeric(R)||~isnumeric(K)
    error('bath_pumped_hahn:generators','Supply conventional or tagged projected generators.');
end
if ~isfield(parameters,'pump')||~isstruct(parameters.pump)||...
   ~all(isfield(parameters.pump,{'center','halfwidth','units','fraction','allow_empty'}))
    error('bath_pumped_hahn:pump','Supply a complete physical pump selector.');
end
end


