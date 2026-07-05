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

% Lorentzian convolution
mex('-R2018a','-O','-DNDBUG','COMPFLAGS=$COMPFLAGS','LINKFLAGS=$LINKFLAGS',...
    [P '/kernel/line_shapes/lorentzcon.cpp'],'-outdir',[P '/kernel/line_shapes']);

% Gaussian convolution
mex('-R2018a','-O','-DNDBUG','COMPFLAGS=$COMPFLAGS','LINKFLAGS=$LINKFLAGS',...
    [P '/kernel/line_shapes/gausscon.cpp'],'-outdir',[P '/kernel/line_shapes']);

% Cubic polynomial roots
mex('-R2018a','-O','-DNDBUG','COMPFLAGS=$COMPFLAGS','LINKFLAGS=$LINKFLAGS',...
    [P '/kernel/eigenfields/cubic_roots.cpp'],'-outdir',[P '/kernel/eigenfields']);

% Sparse matrix assembly
if ispc
    mex('-R2018a','-O','-DNDEBUG','COMPFLAGS=$COMPFLAGS /openmp',...
        [P '/kernel/indexing/fsparse.cpp'],'-outdir',[P '/kernel/indexing']);
elseif isunix&&(~ismac)
    mex('-R2018a','-O','-DNDEBUG','CXXFLAGS=$CXXFLAGS -fopenmp',...
        'LDFLAGS=$LDFLAGS -fopenmp',[P '/kernel/indexing/fsparse.cpp'],...
        '-outdir',[P '/kernel/indexing']);
else
    mex('-R2018a','-O','-DNDEBUG',...
        [P '/kernel/indexing/fsparse.cpp'],'-outdir',[P '/kernel/indexing']);
end

% Sparse double row sorter
mex('-R2018a','-O','-DNDEBUG',...
    [P '/kernel/indexing/spsortrows.cpp'],'-outdir',[P '/kernel/indexing']);

% Sparse double unique columns
mex('-R2018a','-O','-DNDEBUG',...
    [P '/kernel/indexing/spunicols.cpp'],'-outdir',[P '/kernel/indexing']);

end

% Audiophiles don't use their equip-
% ment to listen to your music. Audio-
% philes use your music to listen to
% their equipment.
%
% Alan Parsons
