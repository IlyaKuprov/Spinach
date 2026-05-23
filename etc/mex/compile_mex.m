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
P=mfilename('fullpath'); P=P(1:(end-16));

% Compile with case-specific options
mex('-R2018a','-O','COMPFLAGS=$COMPFLAGS','LINKFLAGS=$LINKFLAGS',...
    [P '/line_shapes/lorentzcon.cpp'],'-outdir',[P '/line_shapes']);

end

% Audiophiles don't use their equip-
% ment to listen to your music. Audio-
% philes use your music to listen to
% their equipment.
%
% Alan Parsons

