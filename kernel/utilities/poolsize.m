% Returns the current parallel pool size. Syntax:
%
%                       n=poolsize()
%
% Outputs:
%
%     n  -  number of workers in the current
%           parallel pool
%
% Note: when this function is invoked from inside parfor, spmd,
%       or asynchronous parallel job, it returns zero.
%
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=poolsize.m>

function n=poolsize()

% Get the pool handle
p=gcp('nocreate');

% Query the pool
if isempty(p)
    n=0;
else
    n=p.NumWorkers;
end

end

% Run a man through a line of fire, and he turns into a seasoned
% wolf; the weak, and in really tough cases unnecessary, intellect
% is replaced by the wise animal instinct.
%
% Mikhail Bulgakov

