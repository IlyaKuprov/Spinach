% Discovers regression tests from their local declarations. Syntax:
%
%                     entries=discover_tests()
%
% Outputs:
%
%     entries - test identifiers, display names, and function names
%               in deterministic identifier order

function entries=discover_tests()

% Search the same test trees that the runner places on the MATLAB path
root=fileparts(fileparts(mfilename('fullpath')));
files=[dir(fullfile(root,'kernel','**','test_*.m'));...
       dir(fullfile(root,'interfaces','**','test_*.m'))];
entries=struct('id',{},'name',{},'function',{});
for n=1:numel(files)

    % Read metadata without executing a test or loading its dependencies
    source=fileread(fullfile(files(n).folder,files(n).name));
    tag=regexp(source,'^% Test: ([^\r\n]+)','tokens','once');
    if isempty(tag)
        error('Expected one Test declaration in %s.',files(n).name);
    end
    if strcmp(tag{1},'manual'), continue; end

    % Keep the public identifier and label independent of the MATLAB filename
    fields=regexp(tag{1},'^(\S+) \| (\S.*)$','tokens','once');
    if isempty(fields)
        error('Malformed Test declaration in %s.',files(n).name);
    end
    [~,func]=fileparts(files(n).name);
    entries(end+1)=struct('id',fields{1},'name',fields{2},'function',func); %#ok<AGROW>

end

% Reject ambiguous dispatch before any tests are run
if numel(unique({files.name}))~=numel(files)||...
   numel(unique({entries.id}))~=numel(entries)
    error('Discovered tests must have unique identifiers and function names.');
end
[~,order]=sort({entries.id});
entries=entries(order);

end


