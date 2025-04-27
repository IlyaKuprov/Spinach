% Gives Matlab a good smack every time MDCS gets its kni-
% ckers in a twist. Syntax:
%
%                        smack()
% 
% This function shuts down the parallel pool, clears the
% workspace, clears the GPUs, and makes sure there are no
% crashed MDCS jobs left over. This function should only
% be used from the command line.
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=smack.m>

function smack()

% Kill the parallel pool
delete(gcp('nocreate'));

% Clear out crashed jobs
myCluster=parcluster('Processes');
delete(myCluster.Jobs);

% Close all handles
fclose('all');

% Clear the workspace
clear('all'); %#ok<CLALL>

% Reset all GPUs
for n=1:gpuDeviceCount
    gpuDevice(n).reset;
end

end

% The best alibi is to be the victim.

