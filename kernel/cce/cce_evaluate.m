% Evaluates one cluster of a cluster correlation expansion: builds the
% retained subsystem and its generators, and calls the pulse sequence.
% Syntax:
%
%   raw=cce_evaluate(model,cluster,configuration,pulse_sequence,parameters)
%
% Parameters:
%
%    model          - CCE model returned by cce_compile
%
%    cluster        - one element of the cluster array
%                     returned by cce_clusters
%
%    configuration  - bath product state record with
%                     bath_ids, bath_rho, and weight
%
%    pulse_sequence - pulse sequence function handle
%
%    parameters     - resolved clexpn parameters
%
% Outputs:
%
%    raw            - record with cluster_id, value, abs_error_bound,
%                     valid, weight, and the diagnostic fields that
%                     the pulse sequence has returned
%
% talos@spindynamics.org

function raw=cce_evaluate(model,cluster,configuration,pulse_sequence,parameters)

% Check consistency
grumble(model,cluster,configuration,pulse_sequence,parameters);

% Build the retained subsystem and its generators
if isfield(configuration,'assignment')
    configuration=rmfield(configuration,'assignment');
end
[spin_system,~]=cce_subsystem(model,parameters,cluster.sites,configuration,[]);
[H,R,K,parameters]=cce_generators(spin_system,parameters,[]);

% Test invariance unless the user declares the diagonal register approximation
if isstruct(H)&&(~isfield(parameters.cce,'projection')||...
                 strcmp(parameters.cce.projection,'invariant'))&&...
   ~parameters.cce.second_order.enabled&&...
   parameters.cce_aux.projection_leakage>spin_system.dimension*eps*H.hamiltonian_norm
    error('cce:population_transfer','Projected states must be invariant under the cluster Hamiltonian.');
end

% Call the pulse sequence
[answer,aux]=pulse_sequence(spin_system,parameters,H,R,K);

% Preserve the natural numeric layout of the pulse sequence output
raw.cluster_id=cluster.id;
raw.value=answer;
raw.abs_error_bound=aux.abs_error_bound;
raw.valid=isfinite(answer);
raw.weight=configuration.weight;
if isfield(parameters.cce_aux,'generator')
    raw.generator=parameters.cce_aux.generator;
end
if isfield(aux,'eigen_residual')
    raw.eigen_residual=aux.eigen_residual;
    raw.orthogonality_residual=aux.orthogonality_residual;
    raw.hamiltonian_bound=aux.hamiltonian_bound;
    raw.backend=aux.backend;
end

% Retain the inclusion decision for each observation time
if isfield(aux,'pulse_time')
    raw.pulse_time=aux.pulse_time;
    raw.pulse_applied=aux.pulse_applied;
end

% Retain generalized echo density
if isfield(aux,'density')
    raw.density=aux.density;
end

% Retain unreduced noise moments and their validity evidence
if isfield(aux,'moment')
    raw.moment=aux.moment;
    raw.mean_time=aux.mean_time;
    raw.mean_zero=aux.mean_zero;
    raw.valid=aux.valid;
end

end

% Consistency enforcement
function grumble(model,cluster,configuration,pulse_sequence,parameters)
if ~isstruct(model)||~isfield(model,'cce')||...
   ~isstruct(cluster)||~isscalar(cluster)||~all(isfield(cluster,{'id','sites'}))||...
   ~isstruct(configuration)||~isscalar(configuration)||~isfield(configuration,'weight')||...
   ~isa(pulse_sequence,'function_handle')||...
   ~isstruct(parameters)||~isscalar(parameters)||~all(isfield(parameters,{'cce','times'}))
    error('cce_evaluate:input','Supply a compiled model, one cluster, one bath configuration, a pulse sequence handle and resolved parameters.');
end
if isfield(parameters.cce,'checkpoint')
    error('cce_options:checkpoint','Only the scheduler may consume checkpoint requests.');
end
dc=prod(model.comp.mults(parameters.observed));
components=1;
if isfield(parameters,'output')&&strcmp(parameters.output,'density')
    components=dc^2;
end
if isfield(parameters,'output')&&strcmp(parameters.output,'observable')&&...
   isfield(parameters,'coil')&&iscell(parameters.coil)
    components=numel(parameters.coil);
end
if 3*numel(parameters.times)*(25*components+16*dc^2)>parameters.cce.numerics.max_memory_bytes
    error('clexpn:memory','Retained time records exceed cce.numerics.max_memory_bytes.');
end
end


