% A shortcut to the parallel pool shutdown command. When
% a Ctrl-C is pressed halfway through a parallel calcula-
% tion, Matlab would sometimes not terminate the worker
% processes, and the parallel pool would get stuck until
% all workers have computed their chunks. Syntax:
%
%                        smack()
% 
% This function forcibly shuts down the parallel pool and
% clears the workspace; do not use it in any way other
% than from the command line.
%
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=smack.m>

function smack()

% Kill the parallel pool
delete(gcp('nocreate'));

% Close all handles
fclose('all');

% Clear variables
clear('all'); %#ok<CLALL>

% Reset GPUs
for n=1:gpuDeviceCount
    gpuDevice(n).reset;
end

end

% The best alibi is to be the victim.

