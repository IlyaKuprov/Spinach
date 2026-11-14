% Read, write or validate a versioned CCE task checkpoint. Syntax:
%
%                  state=cce_checkpoint(filename,action,state)
%
% filename is a MAT file; action is 'read', 'write' or 'validate'. The
% state contains version=1, identity (model/parameter/cluster/schedule/
% sample hashes), payload, task_ids, completed, records and checksum.
% Read requires the expected identity and task_ids in state. Validation
% checks the full schema and checksum without modifying the file. Write
% closes a same-directory temporary before atomic POSIX replacement;
% an existing file with different identity is never overwritten.
%
% Inputs:
%   filename - MAT checkpoint path.
%   action - read, write, or validate.
%   state - version 1 identity, payload, task_ids, completed, records and
%   checksum; for reading, expected identity and task_ids.
%
% Outputs:
%   state - validated/read/written checkpoint state.
%
% See also cce_schedule.
%
% ilya.kuprov@weizmann.ac.il

function state=cce_checkpoint(filename,action,state)

% Validate the operation before touching persistent data
grumble(filename,action,state);
if strcmp(action,'read')

    % Validate stored data and compare all run identities and task IDs
    expected=state; loaded=load(filename,'state');
    if ~isfield(loaded,'state'), error('cce_checkpoint:schema','Missing checkpoint state.'); end
    state=cce_checkpoint(filename,'validate',loaded.state);
    if ~isequaln(state.identity,expected.identity)||~isequaln(state.task_ids,expected.task_ids)
        error('cce_checkpoint:identity','Checkpoint identity or task list differs.');
    end
    return
end

% Require a consistent indexed completion ledger and versioned payload
fields={'version','identity','payload','task_ids','completed','records'};
if ~all(isfield(state,fields))||~isequal(state.version,1)||...
   ~isstruct(state.identity)||~isscalar(state.identity)||...
   ~all(isfield(state.identity,{'model','parameters','clusters','schedule','samples','code'}))||...
   ~isstruct(state.payload)||~isscalar(state.payload)||...
   ~isnumeric(state.task_ids)||size(state.task_ids,2)~=3||...
   any(~isfinite(state.task_ids),'all')||...
   size(unique(state.task_ids,'rows'),1)~=size(state.task_ids,1)||...
   ~islogical(state.completed)||~iscolumn(state.completed)||...
   ~iscell(state.records)||~iscolumn(state.records)||...
   numel(state.completed)~=size(state.task_ids,1)||...
   numel(state.records)~=numel(state.completed)||...
   any(cellfun(@isempty,state.records(state.completed)))
    error('cce_checkpoint:schema','Invalid versioned checkpoint ledger.');
end

% Cover the payload, identities, completion bitmap and raw numerical data
content=state;
if isfield(content,'checksum'), content=rmfield(content,'checksum'); end
checksum=md5_hash(content);
if strcmp(action,'validate')
    if ~isfield(state,'checksum')||~isequal(state.checksum,checksum)
        error('cce_checkpoint:checksum','Checkpoint checksum differs.');
    end
    return
end
state.checksum=checksum;

% Refuse replacement of a checkpoint belonging to another calculation
if isfile(filename)
    existing=cce_checkpoint(filename,'read',state);
    shared=existing.completed&state.completed;
    if ~isequaln(existing.records(shared),state.records(shared))
        error('cce_checkpoint:completed','A completed task cannot be silently overwritten.');
    end
end
folder=fileparts(filename);
if isempty(folder), folder=pwd; end
if ~isfolder(folder), error('cce_checkpoint:directory','Checkpoint directory does not exist.'); end

% Write and verify a closed temporary file before a single atomic rename
temporary=[tempname(folder) '.mat'];
try
    save(temporary,'state','-v7');
    cce_checkpoint(temporary,'read',state);
    [ok,message]=movefile(temporary,filename,'f');
    if ~ok, error('cce_checkpoint:write','%s',message); end
catch exception
    if isfile(temporary), delete(temporary); end
    rethrow(exception);
end

end

% Validate the public call and the expected read identity
function grumble(filename,action,state)
if ~ischar(filename)||~isrow(filename)||isempty(filename)
    error('cce_checkpoint:input','filename must be a nonempty character row.');
end
if ~ischar(action)||~ismember(action,{'read','write','validate'})
    error('cce_checkpoint:input','action must be read, write or validate.');
end
if ~isstruct(state)||~isscalar(state)
    error('cce_checkpoint:input','state must be a scalar structure.');
end
if strcmp(action,'read')&&~all(isfield(state,{'identity','task_ids'}))
    error('cce_checkpoint:input','Read requires expected identity and task_ids.');
end
end

