% Forces a wipe of the Spinach cache folder. Syntax:
% 
%                     wipe_cache(spin_system)
%
% Parameters:
%
%     spin_system - Spinach object with information (stored
%                   in spin_system.sys.scratch) about the
%                   cache folder location, use bootstrap()
%                   to get the default object
%
% Output:
%
%     an attempt is made to delete all all Spinach-specific
%     files in spin_system.sys.scratch; this would fail qui-
%     etly if file system permissions are insufficient
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=wipe_cache.m>

function wipe_cache(spin_system)

% Check consistency
grumble(spin_system);

% Inform the user
report(spin_system,'cache wipe requested by the user...')

% Set cache memory horizon to zero
spin_system.tols.cache_mem=0;

% Call cache management
cacheman(spin_system);

end

% Consistency enforcement
function grumble(spin_system)
if (~isfield(spin_system,'sys'))||...
   (~isfield(spin_system.sys,'scratch'))
    error('scratch folder location information is missing.');
end
if ~exist(spin_system.sys.scratch,'dir')
    error('the specified scratch folder does not exist.');
end
end

% I can't lie to you about your chances, but... you
% have my sympathies.
%
% Bishop the Android,
% in the Alien film

