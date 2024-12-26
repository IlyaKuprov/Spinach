% Returns true if executed inside a parfor or spmd block. This 
% function is used in the internal decision making of Spinach 
% kernel: certain algorithms are switched to their serial ver-
% sions if the calculation is already running inside some par-
% allel loop. Syntax:
%
%                    answer=isworkernode()
%
% Outputs:
%
%    answer - true if running on a parallel worker process
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=isworkernode.m>

function answer=isworkernode()

    % Undocumented function, c'est la vie
    answer=parallel.internal.pool.isPoolWorker();
    
end

% In the beginning the Universe was created. This has 
% made a lot of people very angry and been widely re-
% garded as a bad move.
%
% Douglas Adams, "Hitchhiker's Guide to the Galaxy"

