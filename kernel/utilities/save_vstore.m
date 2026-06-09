% Saves the current parallel pool ValueStore into a Matlab file.
% The snapshot contains keys and values only; callback functions
% are session-local and are not stored. Syntax:
%
%                  save_vstore(file_name)
%
% Parameters:
%
%    file_name - a character string specifying the destination
%                MAT file
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=save_vstore.m>

function save_vstore(file_name)

% Check consistency
grumble(file_name);

% Get the current parallel pool
current_pool=gcp('nocreate');
if isempty(current_pool)
    error('no current parallel pool found.');
end

% Get the current ValueStore
store=current_pool.ValueStore;

% Get all keys and values
key_set=keys(store);
if isempty(key_set)
    val_set=cell(size(key_set));
else
    val_set=get(store,key_set);
end

% Save the snapshot
save(file_name,'key_set','val_set','-v7.3'); drawnow;

end

% Consistency enforcement
function grumble(file_name)
if (~ischar(file_name))||(~isrow(file_name))||isempty(file_name)
    error('file_name must be a non-empty character string.');
end
end

% History has shown us that it's not religion that's
% the problem, but any system of thought that insists
% that one group of people are inviolably in the right,
% whereas the others are in the wrong and must somehow
% be punished.
%
% Rod Liddle

