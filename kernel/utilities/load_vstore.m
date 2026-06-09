% Loads the current parallel pool ValueStore from a Matlab file.
% The current store is cleared before the saved keys and values
% are inserted. Callback functions are session-local and are not
% loaded. Syntax:
%
%                  load_vstore(file_name)
%
% Parameters:
%
%    file_name - a character string specifying the source MAT file
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=load_vstore.m>

function load_vstore(file_name)

% Check consistency
grumble(file_name);

% Load the snapshot
snapshot=load(file_name,'key_set','val_set');

% Check snapshot format
if (~isfield(snapshot,'key_set'))||(~isfield(snapshot,'val_set'))||...
   (~isstring(snapshot.key_set))||(~iscell(snapshot.val_set))||...
   (~isequal(size(snapshot.key_set),size(snapshot.val_set)))
    error('the file must contain a ValueStore snapshot from save_vstore.');
end

% Get the current parallel pool
current_pool=gcp('nocreate');
if isempty(current_pool)
    error('no current parallel pool found.');
end

% Get the current ValueStore
store=current_pool.ValueStore;

% Remove current keys
old_keys=keys(store);
if ~isempty(old_keys)
    remove(store,old_keys);
end

% Insert saved keys and values
if ~isempty(snapshot.key_set)
    put(store,snapshot.key_set,snapshot.val_set);
end

end

% Consistency enforcement
function grumble(file_name)
if (~ischar(file_name))||(~isrow(file_name))||isempty(file_name)
    error('file_name must be a non-empty character string.');
end
if ~isfile(file_name)
    error('file_name does not point to an existing file.');
end
end


