% An include that writes the report of the profiling infrastructure
% around parallel stages. Should be invoked just after a parfor or
% an spmd for which parallel_profiler_start was previously called.
%
% ilya.kuprov@weizmann.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=parallel_profiler_report.m>

% Brief parallel profiler report
if ~isworkernode
    nbytes=mean(tocBytes(gcp),1)/2^20; walltime=toc();
    report(spin_system,['average worker process received ' num2str(nbytes(1)) ...
                        ' MB and sent back ' num2str(nbytes(2)) ' MB']);
    report(spin_system,['parallel stage run time: ' num2str(walltime) ' seconds']);
end

% Detailed parallel profiler report
if (~isworkernode)&&ismember('dafuq',spin_system.sys.enable)
    parpool_history=parProfiler.drainLog(); a=dbstack;
    filename=[spin_system.sys.scratch filesep ...
              datestr(clock,30) '_' a(end-1).name '.mat']; %#ok<CLOCK,DATST> 
    save(filename,'parpool_history');
    report(spin_system,['dafuq saved as ' filename]);
end

% They fuck you up, your mum and dad.   
%     They may not mean to, but they do.   
% They fill you with the faults they had
%     And add some extra, just for you.
% 
% But they were fucked up in their turn
%     By fools in old-style hats and coats,   
% Who half the time were soppy-stern
%     And half at one another's throats.
% 
% Man hands on misery to man.
%     It deepens like a coastal shelf.
% Get out as early as you can,
%     And don't have any kids yourself.
%
% Philip Larkin

% #NHEAD #NGRUM