% Lists Spinach regression tests. Syntax:
%
%                    manifest=list_tests(varargin)
%
% Parameters:
%
%     varargin  - optional name-value pair 'pattern', string
%
% Outputs:
%
%     manifest  - structure array with test identifiers and names
%
% ilya.kuprov@weizmann.ac.il

function manifest=list_tests(varargin)

% Add the test library to the path
root_dir=fileparts(mfilename('fullpath'));
addpath(fullfile(root_dir,'lib'));

% Parse options
options=test_options(varargin{:});
manifest=test_manifest();

% Apply substring filter
if ~isempty(options.pattern)
    keep=contains({manifest.id},options.pattern)|...
         contains({manifest.name},options.pattern);
    manifest=manifest(keep);
end

% Print the list
for n=1:numel(manifest)
    fprintf('%s\t%s\n',manifest(n).id,manifest(n).name);
end

end
