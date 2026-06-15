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
            if islogical(varargin{n+1})&&isscalar(varargin{n+1})
                options.verbose=varargin{n+1};
            elseif ischar(varargin{n+1})&&strcmp(varargin{n+1},'true')
                options.verbose=true;
            elseif ischar(varargin{n+1})&&strcmp(varargin{n+1},'false')
                options.verbose=false;
            elseif isstring(varargin{n+1})&&isscalar(varargin{n+1})&&(varargin{n+1}=="true")
                options.verbose=true;
            elseif isstring(varargin{n+1})&&isscalar(varargin{n+1})&&(varargin{n+1}=="false")
                options.verbose=false;
            else
                error('verbose option must be true or false.');
            end
        case 'stop_on_fail'
            if islogical(varargin{n+1})&&isscalar(varargin{n+1})
                options.stop_on_fail=varargin{n+1};
            elseif ischar(varargin{n+1})&&strcmp(varargin{n+1},'true')
                options.stop_on_fail=true;
            elseif ischar(varargin{n+1})&&strcmp(varargin{n+1},'false')
                options.stop_on_fail=false;
            elseif isstring(varargin{n+1})&&isscalar(varargin{n+1})&&(varargin{n+1}=="true")
                options.stop_on_fail=true;
            elseif isstring(varargin{n+1})&&isscalar(varargin{n+1})&&(varargin{n+1}=="false")
                options.stop_on_fail=false;
            else
                error('stop_on_fail option must be true or false.');
            end
        otherwise
            error(['unknown option: ' varargin{n}]);
    end
end
end
