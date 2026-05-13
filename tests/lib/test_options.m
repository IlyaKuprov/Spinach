% Parses name-value options for the Spinach test runner. Syntax:
%
%                    options=test_options(varargin)
%
% Parameters:
%
%     varargin  - name-value option pairs
%
% Outputs:
%
%     options   - options structure
%
% ilya.kuprov@weizmann.ac.il

function options=test_options(varargin)

% Set defaults
options.pattern='';
options.verbose=false;
options.stop_on_fail=false;

% Parse name-value pairs
if mod(numel(varargin),2)~=0
    error('options must be supplied as name-value pairs.');
end
for n=1:2:numel(varargin)
    switch varargin{n}
        case 'pattern'
            options.pattern=varargin{n+1};
        case 'verbose'
            options.verbose=logical(varargin{n+1});
        case 'stop_on_fail'
            options.stop_on_fail=logical(varargin{n+1});
        otherwise
            error(['unknown option: ' varargin{n}]);
    end
end

end
