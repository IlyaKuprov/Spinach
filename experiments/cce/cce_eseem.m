% Two-pulse ESEEM echo; ideal pi/2 preparation and midpoint pi pulse.
% Syntax:
%
%                 [E,aux]=cce_eseem(spin_system,parameters,H,R,K)
%
% Parameters:
%
%    spin_system - native retained subsystem from cce_subsystem
%
%    parameters  - experiment structure; times is a real vector of
%                  increasing nonnegative total durations in seconds
%
%    H,R,K       - conventional angular-frequency generators, or H
%                  containing explicitly tagged projected branches
%
% transition selects the observed two-state transition. rho0/coil or
% preparation may override the prepared density and detector. normalized
% defaults to true and divides by the detected initial reference.
% Ideal refocusing acts on the selected transition about x. Odd projected
% inversion counts restore the physical coherence before detection.
%
% Outputs:
%
%    E - complex signal array, with one row per requested time and one
%        column per receiver; dimensionless when normalised
%
%    aux.abs_error_bound - time-by-receiver floating-point propagation
%                          estimate, not a rigorous forward error bound
%
% ilya.kuprov@weizmann.ac.il

function [E,aux]=cce_eseem(spin_system,parameters,H,R,K)

% Set the declared signal normalisation default
if ~isfield(parameters,'normalized')
    parameters.normalized=true;
end

% Check consistency
grumble(spin_system,parameters,H,R,K);

% Prepare the density and receiver for the declared transition
[rho0,coil,reference,meta]=cce_observation(spin_system,parameters);
if ~iscell(coil)
    coil={coil};
end

% Allocate the signal and choose the propagation representation
projected=isstruct(H);
E=complex(zeros(numel(parameters.times),numel(reference)));
scale=cellfun(@(x)norm(x,'fro'),coil(:).')./abs(meta.normalization_divisor);
aux.abs_error_bound=zeros(size(E));
aux.error_bound_kind='floating_point_estimate';
if projected
    generator=H;
else
    generator=H+1i*R+1i*K;
end

% Build ideal observed-spin transition pulses without a global bath basis
model=spin_system.model;
obs_dim=spin_system.observed_dimension;
bath_dim=spin_system.bath_dimension;
if ~isfield(model,'transition')
    error('cce_eseem:transition',...
          'An explicit transition is required for ideal refocusing.');
end
alpha=model.transition.alpha;
beta=model.transition.beta;
rotation=eye(obs_dim)-alpha*alpha'-beta*beta'-1i*(alpha*beta'+beta*alpha');
rotation=kron(rotation,eye(bath_dim));
fractions=0.5;

% At zero elapsed time no control is applied, matching the reference signal
for n=1:numel(parameters.times)

    % Start each time point from the same prepared density
    rho=rho0;
    previous=0;
    swapped=false;
    state_error=0;
    if parameters.times(n)>0
        for j=1:numel(fractions)

            % Evolve up to the next ideal inversion
            segment=generator;
            if projected&&swapped
                segment.left=H.right;
                segment.right=H.left;
            end
            duration=(fractions(j)-previous)*parameters.times(n);
            [rho,estimate]=cce_evolution(spin_system,segment,[],rho,duration,1,'final');
            state_error=state_error+estimate;

            % Invert the transition or exchange its projected branches
            if projected
                swapped=~swapped;
            else
                rho=reshape(rho,obs_dim*bath_dim,obs_dim*bath_dim);
                rho=rotation*rho*rotation';
                if strcmp(spin_system.formalism,'zeeman-liouv')
                    rho=rho(:);
                end
            end
            previous=fractions(j);

        end
    end

    % Propagate from the final pulse to the detection time
    segment=generator;
    if projected&&swapped
        segment.left=H.right;
        segment.right=H.left;
    end
    duration=(1-previous)*parameters.times(n);
    [rho,estimate]=cce_evolution(spin_system,segment,[],rho,duration,1,'final');
    aux.abs_error_bound(n,:)=(state_error+estimate)*scale;

    % Restore the detected spin-manifold coherence after an odd inversion count
    if projected&&swapped
        rho=reshape(reshape(rho,bath_dim,bath_dim)',bath_dim^2,1);
    end

    % Detect and normalise each requested receiver
    for k=1:numel(coil)
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

% Validate the registered sequence contract before preparation
function grumble(spin_system,parameters,H,R,K)
if ~isstruct(spin_system)||~isfield(spin_system,'kind')||...
   ~strcmp(spin_system.kind,'cce_subsystem')||~isstruct(parameters)||...
   ~isfield(parameters,'times')||~isnumeric(parameters.times)||...
   ~isreal(parameters.times)||~isvector(parameters.times)||...
   any(~isfinite(parameters.times))||any(parameters.times<0)||any(diff(parameters.times)<0)
    error('cce_eseem:input',...
          'A native retained subsystem and increasing nonnegative times are required.');
end
if isfield(parameters,'sequence')
    error('cce_eseem:sequence',...
          'Use finite_drive, refocusing_train, or this named experiment, not parameters.sequence.');
end
if isfield(parameters,'output')&&strcmp(parameters.output,'density')
    error('cce_eseem:density','This experiment returns detected signal or observables.');
end
if ~(isnumeric(H)||isstruct(H))||~isnumeric(R)||~isnumeric(K)
    error('cce_eseem:generators','Supply conventional or tagged projected generators.');
end
end


