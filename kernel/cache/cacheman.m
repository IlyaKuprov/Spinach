% Cache management heuristics. Looks after the scratch folder and
% prevents it from filling up the disk. Do not call directly.
%
% The function inspects the scratch folder and deletes any files
% that are older than the threshold (default is 365 days) speci-
% fied in spin_system.tols.cache_mem field.
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=cacheman.m>

function cacheman(spin_system)

% Check consistency
grumble(spin_system);

% Calculate the time horizon
time_horizon=now-spin_system.tols.cache_mem; %#ok<TNOW1> 

try % Enforce the ability to write into the scratch directory
    test=1; save([spin_system.sys.scratch filesep 'test.mat'],'test');
    drawnow; delete([spin_system.sys.scratch filesep 'test.mat']); 
catch
    error(['Unable to write into ' spin_system.sys.scratch]);
end

% Look into the scratch directory
dir_cont=dir([spin_system.sys.scratch filesep 'spinach_*']);
if any([dir_cont.datenum]<time_horizon)
    report(spin_system,'cleaning up disk cache, this can take a few seconds...');
end

% Delete anything that is out of date
n_files_gone=0; n_dirs_gone=0;
parfor n=1:numel(dir_cont) % rawrrrrrrr
    if dir_cont(n).datenum<time_horizon
        if dir_cont(n).isdir
            try %#ok<TRYNC> - fail quietly
                rmdir([dir_cont(n).folder filesep dir_cont(n).name],'s');
                n_dirs_gone=n_dirs_gone+1;
            end
        else
            try %#ok<TRYNC> - fail quietly
                delete([dir_cont(n).folder filesep dir_cont(n).name]);
                n_files_gone=n_files_gone+1;
            end
        end
    end
end

% Report to the user
if n_files_gone>0
    report(spin_system,[num2str(n_files_gone) ' out-of-date cache files deleted']);
end
if n_dirs_gone>0
    report(spin_system,[num2str(n_dirs_gone)  ' out-of-date cache directories deleted']);
end

end

% Consistency enforcement
function grumble(spin_system)
if (~isfield(spin_system,'sys'))||(~isfield(spin_system.sys,'scratch'))
    error('the spin_system object does not specify scratch location.');
end
if ~exist(spin_system.sys.scratch,'dir')
    report(spin_system,'expected scratch directory location:');
    report(spin_system,spin_system.sys.scratch);
    error('the scratch directory does not appear to exist.');
end
end

% This model was suggested to Ising by his thesis adviser, Lenz. Ising
% solved the one-dimensional model [...] and on the basis of the fact
% that the one-dimensional model had no phase transition, he asserted
% that there was no phase transition in any dimension. As we shall see,
% this is false. It is ironic that on the basis of an elementary cal-
% culation and erroneous conclusion, Ising's name has become among the
% most commonly mentioned in the theoretical physics literature. But
% history has had its revenge. Ising's name, which is correctly prono-
% unced "Eee-sing", is almost universally mispronounced "Eye-sing".
%
% Barry Simon

