% An include that starts profiling infrastructure around parallel
% stages. Should be invoked just before a parfor or an spmd.
%
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=parallel_profiler_start.m>

% Brief parallel profiler start
if ~isworkernode
    ticBytes(gcp);
end
tic();

% Detailed parallel profiler start
if (~isworkernode)&&ismember('dafuq',spin_system.sys.enable)
    parProfiler=parallel.internal.profiling.PoolProfiler();
end

% In late 1700s, a teacher in a German school asked a kid to
% sum up the numbers from 1 to 100 as a punishment for misbe-
% having. The teacher was astonished when the kid solved the
% problem in seconds:
%
%    S = 1   + 2   + ... + 100
%    S = 100 + 99  + ... + 1
%   --------------------------
%   2S = 101 + 101 + ... + 101   =>  S = 101*100/2 = 5050
%
% The kid's name was Carl Friedrich Gauss.

% #NHEAD #NGRUM