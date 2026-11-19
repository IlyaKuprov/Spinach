% Compile retained-cluster preparation, detection, and normalisation reference.
% Syntax:
%
%   [rho0,coil,reference,meta]=cce_observation(cluster_system,parameters)
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
% detector; generalised output follows cluster_system.formalism.
%
% Parameters:
%
%   cluster_system - retained subsystem returned by cce_subsystem
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

function [rho0,coil,reference,meta]=cce_observation(cluster_system,parameters)

% Validate physical ownership before constructing preparation and detection
grumble(cluster_system,parameters);
model=cluster_system.model;
dc=cluster_system.observed_dimension;
db=cluster_system.bath_dimension;
ids=model.observed_ids;
local=cluster_system;
local.physical_ids=ids;
local.sites=[];
local.multiplicities=model.multiplicities(ids);
local.dimension=dc;
local.isotopes=model.isotopes(ids);
local.physical_to_local=zeros(size(model.physical_ids));
local.physical_to_local(ids)=1:numel(ids);
local.formalism='zeeman-hilb';
local.local_operators=cellfun(@pauli,num2cell(local.multiplicities),...
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
    tolerance=dc*eps;
    if ~isnumeric(observed)||~isequal(size(observed),[dc dc])||...
       any(~isfinite(observed),'all')||norm(observed-observed','fro')>tolerance||...
       abs(trace(observed)-1)>tolerance||min(eig(full((observed+observed')/2)))<-tolerance
        error('cce_observation:density','An explicit preparation is a positive trace-one density matrix.');
    end
elseif isfield(parameters,'preparation')&&isequal(parameters.preparation,'transverse')
    observed=(eye(dc)+2*cce_operator(local,'Lx',ids))/2;
elseif isfield(parameters,'preparation')&&isstruct(parameters.preparation)
    spec=parameters.preparation;
    observed=cce_state(local,spec.states,spec.spins);
elseif isfield(model,'transition')
    vector=(model.transition.alpha+model.transition.beta)/sqrt(2);
    observed=vector*vector';
else
    error('cce_observation:preparation','Declare a preparation or select an observed transition.');
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
elseif isfield(model,'transition')
    detectors=model.transition.alpha*model.transition.beta';
else
    error('cce_observation:detection','Declare a detector or select an observed transition.');
end
if ~iscell(detectors)
    detectors={detectors};
end
if isempty(detectors)||~isvector(detectors)||...
   ~all(cellfun(@(x)isnumeric(x)&&isequal(size(x),[dc dc])&&...
                   all(isfinite(x),'all'),detectors))
    error('cce_observation:detection','Every detector must be a finite observed-spin matrix.');
end
for k=1:numel(detectors)
    detectors{k}=exp(1i*phase)*detectors{k};
end
reference=cellfun(@(x)trace(x'*observed),detectors(:).');
if ~isfield(parameters,'coil')&&~isfield(parameters,'detection')
    reference=exp(-1i*phase)*(model.transition.alpha'*observed*model.transition.beta);
end
if normalised&&any(abs(reference)<=cellfun(@(x)dc*eps*norm(x,'fro')*...
                                             norm(observed,'fro'),detectors(:).'))
    error('cce_observation:reference','Cannot normalize by a vanishing detected reference.');
end

% Bound the requested state and detector allocations before tensor products
if 16*cluster_system.dimension^2*(4+numel(detectors))>...
   model.parameters.cce.numerics.max_memory_bytes
    error('cce:memory','Observation construction exceeds max_memory_bytes.');
end

% Build the bath product in the exact retained order supplied by the subsystem
bath=1;
for k=1:numel(cluster_system.bath_rho)
    bath=kron(bath,cluster_system.bath_rho{k});
end
projected=strcmp(parameters.cce.approximation,'projected');
if projected
    if ~isfield(model,'transition')||numel(detectors)~=1
        error('cce_observation:projected','Projected evolution requires one selected-transition detector.');
    end
    alpha=model.transition.alpha;
    beta=model.transition.beta;
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
    if strcmp(cluster_system.formalism,'zeeman-liouv')
        rho0=rho0(:);
        coil=cellfun(@(x)x(:),coil,'UniformOutput',false);
    end
    if isscalar(coil)
        coil=coil{1};
    end
end

% Keep actual and applied references distinct so normalisation occurs once
meta=struct('physical_ids',cluster_system.physical_ids,'observed_ids',ids,...
            'receiver_phase',phase,'normalized',normalised,...
            'initial_signal',reference,'normalization_divisor',ones(size(reference)),...
            'output_orientation','time-by-detector','vectorization','column-major');
if normalised
    meta.normalization_divisor=reference;
end
if isfield(model,'transition')
    meta.transition=model.transition;
end

end

% Reject ambiguous preparation, detection, and receiver conventions
function grumble(system,parameters)
if isfield(parameters,'sequence')
    error('cce_observation:sequence','Use finite_drive, refocusing_train, or hahn_echo, not parameters.sequence.');
end
if ~isstruct(system)||~isfield(system,'kind')||~strcmp(system.kind,'cce_subsystem')||...
   ~isstruct(parameters)||~isscalar(parameters)||~isfield(parameters,'cce')
    error('cce_observation:input','A native retained subsystem and parameters are required.');
end
if ~isfield(parameters.cce,'approximation')||...
   ~ismember(parameters.cce.approximation,{'generalized','projected'})
    error('cce_observation:approximation','Declare generalized or projected evolution.');
end
if ~ismember(system.formalism,{'zeeman-hilb','zeeman-liouv'})
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
            valid=strcmp(value,'transition')&&isfield(system.model,'transition');
            if strcmp(name{1},'preparation')&&strcmp(value,'transverse')
                valid=isscalar(system.model.observed_ids)&&system.observed_dimension==2;
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
end


