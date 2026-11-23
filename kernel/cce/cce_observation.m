% Compile retained-cluster preparation, detection, and normalisation reference.
% Syntax:
%
%   [rho0,coil,reference,meta]=cce_observation(spin_system,parameters)
%
% Explicit rho0 is an observed-spin density matrix; coil is an observed
% detector matrix or cell vector. Otherwise preparation is 'transition' (the
% equal alpha/beta superposition) or a structure with states and spins for
% cce_state; 'transverse' prepares (identity+2*Lx)/2 for one spin one-half.
% detection is 'transition' (alpha*beta') or a structure with
% operators and spins for cce_operator. Descriptor spins are physical IDs.
% A selected transition supplies the standard preparation/detection only
% when neither an explicit override nor a descriptor is present. No
% polarised state is guessed when no transition or preparation is given.
% receiver_phase is radians and multiplies the detected signal by
% exp(-1i*receiver_phase). normalized is explicit logical; default false.
% reference is the actual initial detected signal before normalisation.
% meta.normalization_divisor is reference if normalized, otherwise ones;
% meta.initial_signal always reports the unnormalised initial detection.
% Nothing in this service divides the state or detector by the reference.
% Projected output is a column-major bath-coherence block with a trace
% detector; generalised output follows spin_system.bas.formalism.
%
% Parameters:
%
%   spin_system - retained subsystem returned by cce_subsystem
%   parameters     - preparation and detection descriptors, optional observed
%                    density rho0 and detector coil, receiver_phase in radians,
%                    and the declared cce.approximation
%
% Outputs:
%
%   rho0      - retained density matrix or column-major density/coherence vector
%   coil      - matching detector matrix, vector, or cell array of detectors
%   reference - initial detected signal, one entry per detector
%   meta      - ordering, receiver phase, and applied normalisation reference
%
% ilya.kuprov@weizmann.ac.il

function [rho0,coil,reference,meta]=cce_observation(spin_system,parameters)

% Validate physical ownership before constructing preparation and detection
grumble(spin_system,parameters);
dc=spin_system.cce.observed_dimension;
db=spin_system.cce.bath_dimension;
ids=parameters.observed;
local=spin_system;
local.cce.ids=ids;
local.cce.sites=[];
local.cce.dimension=dc;
local.cce.physical_to_local=zeros(1,spin_system.comp.nspins);
local.cce.physical_to_local(ids)=1:numel(ids);
local.bas.formalism='zeeman-hilb';
local.cce.local_operators=cellfun(@pauli,num2cell(local.comp.mults(local.cce.ids)),...
                            'UniformOutput',false);
phase=0;
normalised=false;
if isfield(parameters,'receiver_phase')
    phase=parameters.receiver_phase;
end
if isfield(parameters,'normalized')
    normalised=parameters.normalized;
end

% Explicit expert density overrides the named experimental preparation
if isfield(parameters,'rho0')
    observed=parameters.rho0;
elseif isfield(parameters,'preparation')&&isequal(parameters.preparation,'transverse')
    observed=(eye(dc)+2*cce_operator(local,'Lx',ids))/2;
elseif isfield(parameters,'preparation')&&isstruct(parameters.preparation)
    spec=parameters.preparation;
    observed=cce_state(local,spec.states,spec.spins);
else
    vector=(spin_system.cce.transition.alpha+spin_system.cce.transition.beta)/sqrt(2);
    observed=vector*vector';
end
if ~isnumeric(observed)||~isequal(size(observed),[dc dc])||any(~isfinite(observed),'all')
    error('cce_observation:preparation','Preparation must be a finite observed-spin matrix.');
end

% Compile non-Hermitian coils with the same Hilbert-Schmidt convention
if isfield(parameters,'coil')
    detectors=parameters.coil;
elseif isfield(parameters,'detection')&&isstruct(parameters.detection)
    spec=parameters.detection;
    detectors=cce_operator(local,spec.operators,spec.spins);
else
    detectors=spin_system.cce.transition.alpha*spin_system.cce.transition.beta';
end
if ~iscell(detectors)
    detectors={detectors};
end
for k=1:numel(detectors)
    detectors{k}=exp(1i*phase)*detectors{k};
end
reference=cellfun(@(x)trace(x'*observed),detectors(:).');
if ~isfield(parameters,'coil')&&~isfield(parameters,'detection')
    reference=exp(-1i*phase)*(spin_system.cce.transition.alpha'*observed*spin_system.cce.transition.beta);
end
if normalised&&any(abs(reference)<=cellfun(@(x)dc*eps*norm(x,'fro')*...
                                             norm(observed,'fro'),detectors(:).'))
    error('cce_observation:reference','Cannot normalize by a vanishing detected reference.');
end

% Build the bath product in the exact retained order supplied by the subsystem
bath=1;
for k=1:numel(spin_system.cce.bath_rho)
    bath=kron(bath,spin_system.cce.bath_rho{k});
end
projected=strcmp(parameters.cce.approximation,'projected');
if projected
    alpha=spin_system.cce.transition.alpha;
    beta=spin_system.cce.transition.beta;
    basis_detector=alpha*beta';
    coefficient=trace(basis_detector'*detectors{1});
    if norm(detectors{1}-coefficient*basis_detector,'fro')>...
       dc*eps*norm(detectors{1},'fro')
        error('cce_observation:projected','The detector must select the propagated transition block.');
    end
    rho0=(alpha'*observed*beta)*bath;
    rho0=rho0(:);
    coil=coefficient*eye(db);
    coil=coil(:);
else
    rho0=kron(observed,bath);
    coil=cell(size(detectors));
    for k=1:numel(detectors)
        coil{k}=kron(detectors{k},eye(db));
    end
    if strcmp(spin_system.bas.formalism,'zeeman-liouv')
        rho0=rho0(:);
        coil=cellfun(@(x)x(:),coil,'UniformOutput',false);
    end
    if isscalar(coil)
        coil=coil{1};
    end
end

% Keep actual and applied references distinct so normalisation occurs once
meta=struct('physical_ids',spin_system.cce.ids,'observed_ids',ids,...
            'receiver_phase',phase,'normalized',normalised,...
            'initial_signal',reference,'normalization_divisor',ones(size(reference)),...
            'output_orientation','time-by-detector','vectorization','column-major');
if normalised
    meta.normalization_divisor=reference;
end
if isfield(spin_system.cce,'transition')
    meta.transition=spin_system.cce.transition;
end

end

% Reject ambiguous preparation, detection, and receiver conventions
function grumble(system,parameters)
if isfield(parameters,'sequence')
    error('cce_observation:sequence','Use finite_drive, refocusing_train, or hahn_echo, not parameters.sequence.');
end
if ~isstruct(system)||~isfield(system,'cce')||~isfield(system.cce,'sites')||...
   ~isstruct(parameters)||~isscalar(parameters)||~isfield(parameters,'cce')
    error('cce_observation:input','A native retained subsystem and parameters are required.');
end
if ~isfield(parameters.cce,'approximation')||...
   ~ismember(parameters.cce.approximation,{'generalized','projected'})
    error('cce_observation:approximation','Declare generalized or projected evolution.');
end
if ~ismember(system.bas.formalism,{'zeeman-hilb','zeeman-liouv'})
    error('cce_observation:formalism','Only complete Zeeman representations are supported.');
end
if isfield(parameters,'rho0')&&isfield(parameters,'preparation')
    error('cce_observation:preparation','Declare a preparation or rho0, not both.');
end
if isfield(parameters,'coil')&&isfield(parameters,'detection')
    error('cce_observation:detection','Declare a detector descriptor or coil, not both.');
end
for name={'preparation','detection'}
    if isfield(parameters,name{1})
        value=parameters.(name{1});
        if ischar(value)
            valid=strcmp(value,'transition')&&isfield(system.cce,'transition');
            if strcmp(name{1},'preparation')&&strcmp(value,'transverse')
                valid=isscalar(parameters.observed)&&system.cce.observed_dimension==2;
            end
        elseif isstruct(value)&&isscalar(value)
            if strcmp(name{1},'preparation')
                field='states';
            else
                field='operators';
            end
            valid=all(isfield(value,{field,'spins'}));
        else
            valid=false;
        end
        if ~valid
            error('cce_observation:descriptor','Use a supported explicit observation descriptor.');
        end
    end
end
if isfield(parameters,'receiver_phase')&&(~isnumeric(parameters.receiver_phase)||...
   ~isreal(parameters.receiver_phase)||~isscalar(parameters.receiver_phase)||~isfinite(parameters.receiver_phase))
    error('cce_observation:phase','receiver_phase must be a finite real radian scalar.');
end
if isfield(parameters,'normalized')&&(~islogical(parameters.normalized)||~isscalar(parameters.normalized))
    error('cce_observation:normalization','normalized must be an explicit logical scalar.');
end
dc=system.cce.observed_dimension;
if isfield(parameters,'rho0')
    rho0=parameters.rho0;
    if ~isnumeric(rho0)||~isequal(size(rho0),[dc dc])||...
       any(~isfinite(rho0),'all')||norm(rho0-rho0','fro')>dc*eps||...
       abs(trace(rho0)-1)>dc*eps||min(eig(full((rho0+rho0')/2)))<-dc*eps
        error('cce_observation:density','An explicit preparation is a positive trace-one density matrix.');
    end
end
if ~isfield(system.cce,'transition')&&~any(isfield(parameters,{'rho0','preparation'}))
    error('cce_observation:preparation','Declare a preparation or select an observed transition.');
end
if ~isfield(system.cce,'transition')&&~any(isfield(parameters,{'coil','detection'}))
    error('cce_observation:detection','Declare a detector or select an observed transition.');
end
count=1;
if isfield(parameters,'coil')
    coil=parameters.coil;
    if ~iscell(coil)
        coil={coil};
    end
    if isempty(coil)||~isvector(coil)||...
       ~all(cellfun(@(x)isnumeric(x)&&isequal(size(x),[dc dc])&&...
                       all(isfinite(x),'all'),coil))
        error('cce_observation:detection','Every detector must be a finite observed-spin matrix.');
    end
    count=numel(coil);
end
if 16*system.cce.dimension^2*(4+count)>system.tols.cce_mem
    error('cce:memory','Observation construction exceeds sys.tols.cce_mem.');
end
if strcmp(parameters.cce.approximation,'projected')&&...
   (~isfield(system.cce,'transition')||count~=1)
    error('cce_observation:projected','Projected evolution requires one selected-transition detector.');
end
end


