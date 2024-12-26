% Forces a wipe of the Spinach cache folder. Syntax:
% 
%                     wipe_cache(spin_system)
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=wipe_cache.m>

function wipe_cache(spin_system)

% Inform the user
report(spin_system,'cache wipe specifically requested by the user...')

% Set cache memory horizon to zero
spin_system.tols.cache_mem=0;

% Call cache management
cacheman(spin_system);

end

% I can't lie to you about your chances, but... you
% have my sympathies.
%
% Bishop the Android,
% in the Alien film

