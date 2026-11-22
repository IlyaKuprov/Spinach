% Return the explicit retained-cluster experiment contract. Syntax:
%
%                 contract=cce_capabilities(pulse_sequence)
%
% pulse_sequence is a registered five-argument experiment handle. The
% contract declares generator representations, reduction algebra, and the
% evaluator used by the context. Unregistered handles are rejected before
% any basis or generator is constructed; names are not matched heuristically.
% Master-equation task records retain generator evidence; capability versions
% distinguish that raw-record schema when validating a checkpoint.
%
% Parameters:
%
%   pulse_sequence - registered five-argument experiment function handle
%
% Outputs:
%
%   contract - supported representations, outputs, reduction algebra,
%              evaluator handle, and checkpoint capability version
%
% ilya.kuprov@weizmann.ac.il

function contract=cce_capabilities(pulse_sequence)

% Require an actual registered function handle
grumble(pulse_sequence);

% These experiments explicitly support conventional and projected generators
contract.pulse_sequence=pulse_sequence;
contract.formalisms={'zeeman-hilb','zeeman-liouv'};
contract.projected=true;
contract.dissipative=false;
contract.reduction='multiplicative';
contract.density_reduction=false;
contract.outputs={'signal','observable'};
contract.evaluator=@evaluate;
contract.version=2;
if isequal(pulse_sequence,@hahn_echo)||isequal(pulse_sequence,@free_evolution)||...
   isequal(pulse_sequence,@finite_drive)||isequal(pulse_sequence,@refocusing_train)
    contract.version=4;
end
if isequal(pulse_sequence,@hahn_echo)
    contract.version=8;
end
if isequal(pulse_sequence,@bath_pumped_hahn)
    contract.version=5;
end
if isequal(pulse_sequence,@free_evolution)
    contract.version=5;
end
if isequal(pulse_sequence,@free_evolution)||isequal(pulse_sequence,@finite_drive)||...
   isequal(pulse_sequence,@hahn_echo)
    contract.dissipative=true;
end
if isequal(pulse_sequence,@free_evolution)||isequal(pulse_sequence,@finite_drive)||isequal(pulse_sequence,@refocusing_train)||...
   isequal(pulse_sequence,@hahn_echo)
    contract.density_reduction=true;
    contract.outputs{end+1}='density';
end
if isequal(pulse_sequence,@finite_drive)||isequal(pulse_sequence,@refocusing_train)
    contract.projected=false;
end
if isequal(pulse_sequence,@refocusing_train)
    contract.formalisms={'zeeman-hilb'};
end
if isequal(pulse_sequence,@noise_correlation)
    contract.projected=false;
    contract.dissipative=true;
    contract.reduction='additive';
    contract.outputs={'correlation'};
    contract.evaluator=@evaluate_noise;
    contract.version=4;
    contract.auxiliary={'moment','mean_time','mean_zero','valid','abs_error_bound'};
end

% Reject checkpoints predating native backend and physical echo-phase corrections
contract.version=contract.version+2;
if isequal(pulse_sequence,@finite_drive)||isequal(pulse_sequence,@free_evolution)||...
   isequal(pulse_sequence,@cpmg_echo)||isequal(pulse_sequence,@cce_eseem)
    contract.version=contract.version+1;
end

end

% Evaluate a noise request only through its registered experiment contract
function raw=evaluate_noise(model,cluster,configuration,schedule,parameters)

% Preserve the noise experiment and output ownership at the direct entry point
if ~isequal(schedule.pulse_sequence,@noise_correlation)||...
   (isfield(parameters,'output')&&~strcmp(parameters.output,'correlation'))
    error('clexpn:noise','Noise evaluation requires noise_correlation and its correlation output.');
end
raw=evaluate(model,cluster,configuration,schedule,parameters);

end

% Evaluate one native subsystem through the ordinary pulse-sequence signature
function raw=evaluate(model,cluster,configuration,schedule,parameters)

% Only the scheduler may consume checkpoint requests
if isfield(parameters.cce,'checkpoint')
    error('cce_options:checkpoint','Direct experiment evaluators cannot own a checkpoint ledger.');
end

% Bound retained time records before allocating subsystem or propagation arrays
dc=prod(model.multiplicities(model.observed_ids));
components=1;
if isfield(parameters,'output')&&strcmp(parameters.output,'density')
    components=dc^2;
end
if isfield(parameters,'output')&&strcmp(parameters.output,'observable')&&...
   isfield(parameters,'coil')&&iscell(parameters.coil)
    components=numel(parameters.coil);
end
bytes=3*numel(parameters.times)*(25*components+16*dc^2);
if ~isfinite(bytes)||bytes>parameters.cce.numerics.max_memory_bytes
    error('clexpn:memory','Retained time records exceed cce.numerics.max_memory_bytes.');
end

% Build only the retained subsystem and its representation-compatible generators
if isequal(schedule.pulse_sequence,@noise_correlation)
    parameters.cce.approximation='generalized';
end
physical=configuration;
if isfield(physical,'assignment')
    physical=rmfield(physical,'assignment');
end
[spin_system,~]=cce_subsystem(model,cluster.sites,physical,[]);
[H,R,K,parameters]=cce_generators(spin_system,parameters,[]);

% Test invariance unless the user declares the diagonal register approximation
if isstruct(H)&&(~isfield(parameters.cce,'projection')||...
                 strcmp(parameters.cce.projection,'invariant'))&&...
   ~parameters.cce.second_order.enabled&&...
   parameters.cce_aux.projection_leakage>spin_system.dimension*eps*H.hamiltonian_norm
    error('cce:population_transfer','Projected states must be invariant under the cluster Hamiltonian.');
end
pulse_sequence=schedule.pulse_sequence;
[answer,aux]=pulse_sequence(spin_system,parameters,H,R,K);

% Preserve the experiment's natural numeric layout in the scheduler record
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

% Retain the experiment-owned inclusion decision for each observation time
if isfield(aux,'pulse_time')
    raw.pulse_time=aux.pulse_time;
    raw.pulse_applied=aux.pulse_applied;
end

% Retain generalized echo density and numerical propagation diagnostics
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

% Validate the registry key without inspecting its textual spelling
function grumble(pulse_sequence)
if ~isa(pulse_sequence,'function_handle')
    error('cce_capabilities:handle','pulse_sequence must be a function handle.');
end
registered={@free_evolution,@finite_drive,@refocusing_train,@hahn_echo,...
            @cpmg_echo,@bath_pumped_hahn,@cce_eseem,@noise_correlation};
if ~any(cellfun(@(sequence)isequal(sequence,pulse_sequence),registered))
    error('cce_capabilities:unregistered','Experiment is not registered for clexpn.');
end
end


