% Bind a context and experiment without constructing full-system generators.
% Syntax:
%
%       task=context_task(context,pulse_sequence,parameters,assumptions)
%
% Parameters:
%
%   context        - @clexpn; the nested context owns its generators
%   pulse_sequence - registered experiment function handle
%   parameters     - scalar experiment parameter structure
%   assumptions    - Hamiltonian assumptions, a character row
%
% Outputs:
%
%   task - deferred context descriptor for orientation_average or powder;
%          the four input fields are stored without evaluating the context
%
% ilya.kuprov@weizmann.ac.il

function task=context_task(context,pulse_sequence,parameters,assumptions)

% Validate the deferred calculation before storing its inputs
grumble(context,pulse_sequence,parameters,assumptions);
task=struct('context',context,'pulse_sequence',pulse_sequence,...
            'parameters',parameters,'assumptions',assumptions);

end

% Validate the supported context and its unevaluated call arguments
function grumble(context,pulse_sequence,parameters,assumptions)
if ~isa(context,'function_handle')||~isequal(context,@clexpn)
    error('context_task:context','Only the generator-owning clexpn context is supported.');
end
if ~isa(pulse_sequence,'function_handle')
    error('context_task:sequence','Supply an experiment function handle.');
end
if ~isstruct(parameters)||~isscalar(parameters)
    error('context_task:parameters','Supply a scalar parameter structure.');
end
if ~ischar(assumptions)||~isrow(assumptions)||isempty(assumptions)
    error('context_task:assumptions','Supply a nonempty assumptions character row.');
end
cce_capabilities(pulse_sequence);
end


