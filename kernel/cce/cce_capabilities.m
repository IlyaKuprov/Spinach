% Return the explicit retained-cluster experiment contract. Syntax:
%
%                 contract=cce_capabilities(pulse_sequence)
%
% pulse_sequence is a registered five-argument experiment handle. The
% contract declares generator representations, reduction algebra and the
% evaluator used by the context. Unregistered handles are rejected before
% any basis or generator is constructed; names are not matched heuristically.
%
% ilya.kuprov@weizmann.ac.il

function contract=cce_capabilities(pulse_sequence)

% Require an actual registered function handle
grumble(pulse_sequence);
registered={@free_evolution,@hahn_echo,@cpmg_echo,@bath_pumped_hahn,@cce_eseem,@noise_correlation};
if ~any(cellfun(@(sequence)isequal(sequence,pulse_sequence),registered))
    error('cce_capabilities:unregistered','Experiment is not registered for clexpn.');
end

% These experiments explicitly support conventional and projected generators
contract.pulse_sequence=pulse_sequence;
contract.formalisms={'zeeman-hilb','zeeman-liouv'};
contract.projected=true; contract.dissipative=false;
contract.reduction='multiplicative'; contract.density_reduction=false;
contract.outputs={'signal','observable'};
contract.evaluator=@evaluate; contract.version=2;
if isequal(pulse_sequence,@free_evolution)||isequal(pulse_sequence,@hahn_echo)
    contract.dissipative=true;
end
if isequal(pulse_sequence,@free_evolution)||isequal(pulse_sequence,@hahn_echo)
    contract.density_reduction=true; contract.outputs{end+1}='density';
end
if isequal(pulse_sequence,@noise_correlation)
    contract.projected=false; contract.dissipative=true;
    contract.reduction='additive'; contract.outputs={'correlation'};
    contract.auxiliary={'moment','mean_time','mean_zero','valid','abs_error_bound'};
end

end

% Evaluate one native subsystem through the ordinary pulse-sequence signature
function raw=evaluate(model,cluster,configuration,schedule,parameters)

% Build only the retained subsystem and its representation-compatible generators
if isequal(schedule.pulse_sequence,@noise_correlation)
    parameters.cce.approximation='generalized';
end
physical=configuration;
if isfield(physical,'assignment'), physical=rmfield(physical,'assignment'); end
[spin_system,~]=cce_subsystem(model,cluster.sites,physical,[]);
[H,R,K,parameters]=cce_generators(spin_system,parameters,[]);
pulse_sequence=schedule.pulse_sequence;
if isequal(pulse_sequence,@noise_correlation)||isequal(pulse_sequence,@free_evolution)
    [answer,aux]=pulse_sequence(spin_system,parameters,H,R,K);
else
    answer=pulse_sequence(spin_system,parameters,H,R,K);
end

% Preserve the experiment's natural numeric layout in the scheduler record
raw.cluster_id=cluster.id; raw.value=answer;
raw.abs_error_bound=zeros(size(answer)); raw.valid=isfinite(answer);
raw.weight=configuration.weight;
if exist('aux','var')&&isfield(aux,'abs_error_bound')
    raw.abs_error_bound=aux.abs_error_bound;
end
if exist('aux','var')&&isfield(aux,'eigen_residual')
    raw.eigen_residual=aux.eigen_residual;
    raw.orthogonality_residual=aux.orthogonality_residual;
    raw.hamiltonian_bound=aux.hamiltonian_bound; raw.backend=aux.backend;
end
if exist('aux','var')&&isfield(aux,'moment')
    raw.moment=aux.moment; raw.mean_time=aux.mean_time; raw.mean_zero=aux.mean_zero;
    raw.abs_error_bound=aux.abs_error_bound; raw.valid=aux.valid;
end

end

% Validate the registry key without inspecting its textual spelling
function grumble(pulse_sequence)
if ~isa(pulse_sequence,'function_handle')
    error('cce_capabilities:handle','pulse_sequence must be a function handle.');
end
end


