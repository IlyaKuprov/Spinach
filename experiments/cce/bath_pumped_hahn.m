% Bath-pumped Hahn echo with an explicitly selected spectral window.
%
% Syntax:
%
%                 [E,aux]=bath_pumped_hahn(spin_system,parameters,H,R,K)
%
% Parameters:
%
%     spin_system - retained subsystem supplied by clexpn
%
%     parameters  - experiment settings described below
%
%     H,R,K       - conventional generators, or the projected branch pair in H
%
% Outputs:
%
%     E   - complex time-by-detector response
%
%     aux - retained propagation and numerical-error diagnostics
%
% parameters.times contains total sequence durations in seconds. transition
% selects the observed two-state transition; rho0/coil or preparation may
% override the prepared density and detector. normalized defaults to true
% and divides by the actual detected initial reference. The pump requires
% pump.center, halfwidth, units (Hz or rad/s), fraction, and allow_empty;
% optional sites/species restrict the bath.
% Ideal refocusing acts on the selected transition about x. H,R,K are the
% conventional generators, or the projected branch pair with its dissipator.
% Detection uses a fixed physical receiver: the observed inversion takes the
% bath-coherence adjoint, including the prepared coherence, before detection.
%
% Pump center and halfwidth are finite scalars, halfwidth is nonnegative;
% fraction is in [0,1] and allow_empty is a logical scalar. sites lists
% distinct physical bath IDs; species is a row cell of isotope labels.
% Exterior selected spins are omitted when cce.meanfield is none: each
% retained cluster evolves its marginal bath state. Static conditioned
% exterior spins cannot be pumped without time-dependent mean fields.
% pump.boundary may be open (strict, legacy physical window) or closed
% (inclusive, the native default). A zero-width open window selects nothing.
% aux.selected_sites records the full physical selection before cluster pruning;
% aux.boundary_rule records the interval convention. Spectral/error diagnostics
% are also retained for coherent projected echoes. On the interval path,
% aux.abs_error_bound estimates free propagation and numerical pulse actions,
% not preparation or readout rounding. For GKSL evolution, the remaining
% channels contract trace norm, so converting each local Frobenius error to
% trace norm costs at most sqrt(d), also for an embedded d-by-d coherence.
%
% ilya.kuprov@weizmann.ac.il

function [E,aux]=bath_pumped_hahn(spin_system,parameters,H,R,K)

% Set the declared signal normalisation default
if ~isfield(parameters,'normalized')
    parameters.normalized=true;
end

% Validate the experiment inputs
grumble(spin_system,parameters,H,R,K);

% Prepare the physical density and receiver
aux=struct();
[rho0,coil,reference,meta]=cce_observation(spin_system,parameters);
if ~iscell(coil)
    coil={coil};
end
projected=iscell(H);
E=complex(zeros(numel(parameters.times),numel(reference)));
if projected
    generator=[H {R}];
else
    generator=H+1i*R+1i*K;
end

% Build ideal observed-spin transition pulses without a global bath basis
dc=spin_system.cce.observed_dimension;
db=spin_system.cce.bath_dimension;
a=spin_system.cce.transition.alpha;
b=spin_system.cce.transition.beta;
rotation=eye(dc)-a*a'-b*b'-1i*(a*b'+b*a');
rotation=kron(rotation,eye(db));

% Select against canonical laboratory-frame bath frequencies in rad/s
spec=parameters.pump;
ids=parameters.cce.bath_spins;
if isfield(spec,'sites')
    ids=intersect(ids,spec.sites,'stable');
end
if isfield(spec,'species')
    ids=ids(ismember(spin_system.comp.isotopes(ids),cellstr(spec.species)));
end
frequencies=zeros(size(ids));
for k=1:numel(ids)

    % Read the axial bath Zeeman frequency
    tensor=spin_system.cce.zeeman{ids(k)};
    frequencies(k)=tensor(3,3);

end
factor=1;
if strcmp(spec.units,'Hz')
    factor=2*pi;
end
distance=abs(frequencies-factor*spec.center);
width=factor*spec.halfwidth;
if isfield(spec,'boundary')&&strcmp(spec.boundary,'open')
    ids=ids(distance<width);
    rule='strict_abs_less_than';
else
    ids=ids(distance<=width);
    rule='abs_less_than_or_equal';
end
if isempty(ids)&&~spec.allow_empty
    error('bath_pumped_hahn:empty','Pump window selected no bath spins.');
end
if isfield(parameters.cce,'meanfield')&&strcmp(parameters.cce.meanfield,'static')&&...
   ~all(ismember(ids,spin_system.cce.sites))
    error('bath_pumped_hahn:exterior','Pumping static exterior mean fields requires time-dependent conditioning.');
end
aux.selected_sites=ids;
aux.boundary_rule=rule;
ids=intersect(ids,spin_system.cce.sites,'stable');
[fractions,order]=sort([0.5 spec.fraction]);
pulse_kind=[1 2];
pulse_kind=pulse_kind(order);

% Retain ordered branch propagators and their numerical error diagnostics
if projected&&isempty(R)&&...
   ~isfield(parameters,'rho0')&&~isfield(parameters,'preparation')&&~isfield(parameters,'coil')
    selected=aux.selected_sites;

    % Resolve both spectra once for the complete observation-time axis
    [va,ea]=eig(full((H{1}+H{1}')/2),'vector');
    [vb,eb]=eig(full((H{2}+H{2}')/2),'vector');
    ra=norm(H{1}*va-va.*ea.','fro');
    rb=norm(H{2}*vb-vb.*eb.','fro');
    oa=norm(va'*va-eye(size(H{1})),'fro');
    ob=norm(vb'*vb-eye(size(H{2})),'fro');
    residual=max(ra,rb);
    orthogonality=max(oa,ob);
    db=spin_system.cce.bath_dimension;

    % Construct the physical ideal rotations in retained spin order
    transition=spin_system.cce.transition;
    pa=kron(transition.alpha,eye(db));
    pb=kron(transition.beta,eye(db));
    observed=kron((transition.alpha*transition.beta'+...
                   transition.beta*transition.alpha')/2,eye(db));
    pumped=zeros(spin_system.cce.dimension);
    local=spin_system;
    local.bas.formalism='zeeman-hilb';
    for id=reshape(ids,1,[])
        pumped=pumped+cce_operator(local,'Lx',id);
    end
    rotations={expm(-1i*(pi*full(observed))),expm(-1i*(pi*full(pumped)))};
    pulses={{pa'*rotations{1}*pb,pb'*rotations{1}*pa},...
            {pa'*rotations{2}*pa,pb'*rotations{2}*pb}};

    % Assemble the retained product bath density
    rho_b=1;
    for id=spin_system.cce.sites

        % Accumulate the retained physical-spin contribution
        where=find(spin_system.cce.configuration.bath_ids==id,1);
        rho_b=kron(rho_b,spin_system.cce.configuration.bath_rho{where});

    end
    E=complex(zeros(numel(parameters.times),1));
    bounds=zeros(size(E));

    % Match the ordered branch products and their arithmetic error estimate
    for n=1:numel(parameters.times)

        % Propagate and detect the current observation time
        time=parameters.times(n);
        ua=eye(db);
        ub=eye(db);
        swapped=false;
        previous=0;
        if time~=0
            for j=1:numel(fractions)

                % Apply the next free interval and physical rotation
                duration=(fractions(j)-previous)*time;
                left=(va.*exp(-1i*duration*ea).')*va';
                right=(vb.*exp(-1i*duration*eb).')*vb';
                if swapped
                    ua=right*ua;
                    ub=left*ub;
                else
                    ua=left*ua;
                    ub=right*ub;
                end
                pair=pulses{pulse_kind(j)};
                if pulse_kind(j)==1
                    if ~swapped
                        ua=pair{2}*ua;
                        ub=pair{1}*ub;
                    else
                        ua=pair{1}*ua;
                        ub=pair{2}*ub;
                    end
                    swapped=~swapped;
                elseif ~swapped
                    ua=pair{1}*ua;
                    ub=pair{2}*ub;
                else
                    ua=pair{2}*ua;
                    ub=pair{1}*ub;
                end
                previous=fractions(j);

            end
            duration=(1-previous)*time;
            left=(va.*exp(-1i*duration*ea).')*va';
            right=(vb.*exp(-1i*duration*eb).')*vb';
            if swapped
                ua=right*ua;
                ub=left*ub;
            else
                ua=left*ua;
                ub=right*ub;
            end
            E(n)=conj(trace(ua*rho_b*ub'));
        else
            E(n)=1;
        end
        operations=28*db^3;
        gamma=operations*eps/(1-operations*eps);
        bounds(n)=2*time*(residual+spin_system.cce.roundoff_bound)+12*orthogonality+...
                  gamma*(1+time*spin_system.cce.hamiltonian_norm);

    end
    if ~meta.normalized
        E=E*reference;
        bounds=bounds*abs(reference);
    end

    % Retain numerical propagation diagnostics
    aux=struct('abs_error_bound',bounds,'eigen_residual',residual,...
               'orthogonality_residual',orthogonality,...
               'hamiltonian_bound',spin_system.cce.roundoff_bound,...
               'backend',spin_system.cce.backend);

    aux.selected_sites=selected;
    aux.boundary_rule=rule;
    return
end

% Embed retained bath rotations for density and coherence-block propagation
pump_operator=zeros(dc*db);
local=spin_system;
local.bas.formalism='zeeman-hilb';
for id=reshape(ids,1,[])

    % Accumulate the retained physical-spin contribution
    pump_operator=pump_operator+cce_operator(local,'Lx',id);

end
if projected
    pump_operator=pump_operator(1:db,1:db);
end
pump_rotation=expm(-1i*pi*pump_operator);

% Convert Frobenius state errors to each detected component
scale=cellfun(@(x)norm(x,'fro'),coil(:).')./abs(meta.normalization_divisor);
if projected
    pulse_dim=db;
else
    pulse_dim=dc*db;
end

% Bound the Frobenius gain of all subsequent trace-norm-contractive channels
if (projected&&nnz(R)>0)||(~projected&&(nnz(R)>0||nnz(K)>0))
    scale=sqrt(pulse_dim)*scale;
end
aux.abs_error_bound=zeros(size(E));
aux.error_bound_kind='floating_point_estimate';

% At zero elapsed time no control is applied, matching the reference signal
for n=1:numel(parameters.times)

    % Propagate and detect the current observation time
    rho=rho0;
    previous=0;
    swapped=false;
    state_error=0;
    if parameters.times(n)>0
        for j=1:numel(fractions)

            % Apply the next free interval and physical rotation
            segment=generator;
            if projected&&swapped
                segment=generator([2 1 3]);
            end
            duration=(fractions(j)-previous)*parameters.times(n);
            [rho,estimate]=cce_evolution(spin_system,segment,[],rho,duration,1,'final');
            state_error=state_error+estimate;

            % Charge the two-sided pulse product, but not a projected branch relabelling
            if (pulse_kind(j)==1&&~projected)||(pulse_kind(j)==2&&~isempty(ids))
                state_error=state_error+16*pulse_dim^3*eps*norm(rho,'fro');
            end
            if pulse_kind(j)==1
                if projected
                    swapped=~swapped;
                else
                    rho=reshape(rho,dc*db,dc*db);
                    rho=rotation*rho*rotation';
                    if strcmp(spin_system.bas.formalism,'zeeman-liouv')
                        rho=rho(:);
                    end
                end
            else
                if projected
                    rho=reshape(rho,db,db);
                    rho=pump_rotation*rho*pump_rotation';
                    rho=rho(:);
                else
                    rho=reshape(rho,dc*db,dc*db);
                    rho=pump_rotation*rho*pump_rotation';
                    if strcmp(spin_system.bas.formalism,'zeeman-liouv')
                        rho=rho(:);
                    end
                end
            end
            previous=fractions(j);

        end
    end
    segment=generator;
    if projected&&swapped
        segment=generator([2 1 3]);
    end
    duration=(1-previous)*parameters.times(n);
    [rho,estimate]=cce_evolution(spin_system,segment,[],rho,duration,1,'final');
    aux.abs_error_bound(n,:)=(state_error+estimate)*scale;

    % Restore the detected spin-manifold coherence after an odd inversion count
    if projected&&swapped
        rho=reshape(reshape(rho,db,db)',db^2,1);
    end
    for k=1:numel(coil)

        % Evaluate the fixed physical receiver
        if isvector(rho)
            value=coil{k}'*rho;
        else
            value=trace(coil{k}'*rho);
        end
        if meta.normalized
            value=value/reference(k);
        end
        E(n,k)=value;

    end

end

end

% Consistency enforcement
function grumble(spin_system,parameters,H,R,K)
if ~isstruct(spin_system)||~isfield(spin_system,'cce')||...
   ~isfield(spin_system.cce,'sites')||~isstruct(parameters)||...
   ~isfield(parameters,'times')||~isnumeric(parameters.times)||...
   ~isreal(parameters.times)||~isvector(parameters.times)||...
   any(~isfinite(parameters.times))||any(parameters.times<0)||any(diff(parameters.times)<0)
    error('bath_pumped_hahn:input','A native retained subsystem and increasing nonnegative times are required.');
end
if isfield(parameters,'sequence')
    error('bath_pumped_hahn:sequence','Use finite_drive, refocusing_train, or this named experiment, not parameters.sequence.');
end
if isfield(parameters,'output')&&strcmp(parameters.output,'density')
    error('bath_pumped_hahn:density','This experiment returns detected signal or observables.');
end
if isfield(parameters,'bath_jumps')&&~isempty(parameters.bath_jumps)
    error('bath_pumped_hahn:dissipation','This experiment is closed-system only.');
end
if ~(isnumeric(H)||iscell(H))||~isnumeric(R)||~isnumeric(K)
    error('bath_pumped_hahn:generators','Supply conventional or tagged projected generators.');
end
if ~isfield(parameters,'pump')||~isstruct(parameters.pump)||...
   ~all(isfield(parameters.pump,{'center','halfwidth','units','fraction','allow_empty'}))
    error('bath_pumped_hahn:pump','Supply a complete physical pump selector.');
end
spec=parameters.pump;
if isfield(spec,'boundary')&&(~ischar(spec.boundary)||~ismember(spec.boundary,{'open','closed'}))
    error('bath_pumped_hahn:pump','The pump boundary must be open or closed.');
end
if ~isscalar(spec)||~isnumeric(spec.center)||~isreal(spec.center)||...
   ~isscalar(spec.center)||~isfinite(spec.center)||...
   ~isnumeric(spec.halfwidth)||~isreal(spec.halfwidth)||...
   ~isscalar(spec.halfwidth)||~isfinite(spec.halfwidth)||spec.halfwidth<0||...
   ~ischar(spec.units)||~ismember(spec.units,{'Hz','rad/s'})||...
   ~isnumeric(spec.fraction)||~isreal(spec.fraction)||...
   ~isscalar(spec.fraction)||~isfinite(spec.fraction)||...
   spec.fraction<0||spec.fraction>1||...
   ~islogical(spec.allow_empty)||~isscalar(spec.allow_empty)
    error('bath_pumped_hahn:pump','Use finite physical units, nonnegative width and a fraction in [0,1].');
end
if isfield(spec,'sites')&&(~isnumeric(spec.sites)||~isreal(spec.sites)||...
   (~isempty(spec.sites)&&~isrow(spec.sites))||...
   any(~ismember(spec.sites,parameters.cce.bath_spins))||...
   numel(unique(spec.sites))~=numel(spec.sites))
    error('bath_pumped_hahn:sites','Pump sites must be distinct physical bath IDs.');
end
if isfield(spec,'species')&&(~iscellstr(spec.species)||~isrow(spec.species)||...
   any(~ismember(spec.species,spin_system.comp.isotopes)))
    error('bath_pumped_hahn:species','Pump species must be a row of isotope labels from the spin_system.');
end
if ~isfield(spin_system.cce,'transition')
    error('bath_pumped_hahn:transition','An explicit transition is required for ideal refocusing.');
end
ids=parameters.cce.bath_spins;
if isfield(spec,'sites')
    ids=intersect(ids,spec.sites,'stable');
end
if isfield(spec,'species')
    ids=ids(ismember(spin_system.comp.isotopes(ids),spec.species));
end
if any(cellfun(@(z)norm(z(3,1:2))>3*eps*norm(z(3,:)),spin_system.cce.zeeman(ids)))
    error('bath_pumped_hahn:axial','Spectral selection requires an axial canonical Zeeman field.');
end
end


