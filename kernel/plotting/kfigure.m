% Resets the stupid ass figure defaults in R2025a 
% and later back to sensible values. Syntax:
%
%                   handle=kfigure(varargin)
%
% Parameters:
%
%    varargin - same arguments as those accepted by
%               Matlab's figure function
%
% Outputs:
%
%    handle   - figure handle
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=kfigure.m>

function handle=kfigure(varargin)

% Check consistency
grumble(varargin{:});

% Reset to pre-R2025a settings
set(groot,'defaultFigurePosition',[680 458 560 420]); 
set(groot,'defaultFigureWindowStyle','normal'); 
set(groot,'defaultFigureMenuBar','figure'); 
set(groot,'defaultFigureToolbar','figure'); 

% Create and return a handle
handle=figure(varargin{:});

end

% Consistency enforcement
function grumble(varargin)
if isempty(varargin)
    return
end
if isnumeric(varargin{1})
    if (~isscalar(varargin{1}))||(~isreal(varargin{1}))||(~isfinite(varargin{1}))
        error('the figure number must be a finite real scalar.');
    end
    varargin=varargin(2:end);
end
if mod(numel(varargin),2)~=0
    error('property/value arguments must come in pairs.');
end
for n=1:2:numel(varargin)
    if (~ischar(varargin{n}))&&(~(isstring(varargin{n})&&isscalar(varargin{n})))
        error('property names must be character strings.');
    end
end
end

% The most common error of a smart engineer is to 
% optimize a thing that should not exist.
%
% Elon Musk


