% Spinach setup script. Nobody ever reads the documentation,
% so hopefully they would see this function and run it. If
% no input arguments are supplied, that means the user did
% not even read this header - oy vey, then we assume a PhD
% student with a laptop. Otherwise, there are a few specific
% config options for different system types. Syntax:
%
%                 fix_path(config_style)
%
% ilya.kuprov@weizmann.ac.il
% 
% <https://spindynamics.org/wiki/index.php?title=fix_path.m>

function fix_path(config_style)

% Default config style
if ~exist('config_style','var')
    config_style='noob';
end

% Check consistency
grumble(config_style);

% Resolve the Spinach root directory
spinach_root=fileparts(mfilename('fullpath'));

% Run the configuration
switch config_style

    case {'noob','reset'}

        % Status report
        disp('Resetting Matlab path...');

        % Reset Matlab path to default
        restoredefaultpath();

        % Status report
        disp('Updating Matlab path...');

        % Add Spinach directories to path
        addpath(genpath(fullfile(spinach_root,'etc')),'-begin');
        addpath(genpath(fullfile(spinach_root,'experiments')),'-begin');
        addpath(genpath(fullfile(spinach_root,'interfaces')),'-begin');
        addpath(genpath(fullfile(spinach_root,'kernel')),'-begin');

        % Run existential checks
        existentials();

        % Report to the console
        disp('Spinach is ready to run.');
       
    case 'add'

        % Status report
        disp('Updating Matlab path...');

        % Keep revious path, add Spinach
        addpath(genpath(fullfile(spinach_root,'etc')),'-begin');
        addpath(genpath(fullfile(spinach_root,'experiments')),'-begin');
        addpath(genpath(fullfile(spinach_root,'interfaces')),'-begin');
        addpath(genpath(fullfile(spinach_root,'kernel')),'-begin');

        % Run existential checks
        existentials();

        % Report to the console
        disp('Spinach is ready to run.');

    case 'remove'

        % Status report
        disp('Updating Matlab path...');

        % Keep revious path, remove Spinach
        rmpath(genpath(fullfile(spinach_root,'etc')));
        rmpath(genpath(fullfile(spinach_root,'experiments')));
        rmpath(genpath(fullfile(spinach_root,'interfaces')));
        rmpath(genpath(fullfile(spinach_root,'kernel')));

        % Report to the console
        disp('Spinach folders have been removed from Matlab path.');

    otherwise

        % Complain and bomb out
        error('unknown configuration style.');

end

end

% Consistency enforcement
function grumble(config_style)
if ~ischar(config_style)
    error('config_style must be a character string.');
end
end

% Listening to Dirac was a dreadful experience. We accepted Dirac's ideas
% and fell under their influence only when we read his papers. But at a
% seminar... Dirac would come in, no smile, no enthusiasm. He would take a
% piece of chalk in his long fingers and start writing formulae on the
% board, without saying a word. After a while, Born couldn't stand it:
% "Paul, tell us, what are you writing?" Dirac, still writing, reluctantly
% started to speak: "W minus alpha ar pi ar minus alpha zero m c and all
% this multiplied by psi, then alpha mu by alpha nu..." and so on in the
% same vein - and he sincerely believed that he was explaining.
%
% Yuri B. ("Georg") Rumer

