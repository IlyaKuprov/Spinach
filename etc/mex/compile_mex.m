% MEX compilation utility. Rebuilds all C++ MEX binaries in 
% the current directory. Syntax:
%
%                       compile_mex()
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=compile_mex.m>

function compile_mex() % #NGRUM #NHEAD

% Get own location
P=mfilename('fullpath'); P=P(1:(end-20));

% Compile with case-specific options
mex('-R2018a','-O','-DNDBUG','COMPFLAGS=$COMPFLAGS','LINKFLAGS=$LINKFLAGS',...
    [P '/kernel/line_shapes/lorentzcon.cpp'],'-outdir',[P '/kernel/line_shapes']);
mex('-R2018a','-O','-DNDBUG','COMPFLAGS=$COMPFLAGS','LINKFLAGS=$LINKFLAGS',...
    [P '/kernel/line_shapes/gausscon.cpp'],'-outdir',[P '/kernel/line_shapes']);
mex('-R2018a','-O','-DNDBUG','COMPFLAGS=$COMPFLAGS','LINKFLAGS=$LINKFLAGS',...
    [P '/kernel/eigenfields/cubic_roots.cpp'],'-outdir',[P '/kernel/eigenfields']);

end

% Audiophiles don't use their equip-
% ment to listen to your music. Audio-
% philes use your music to listen to
% their equipment.
%
% Alan Parsons

