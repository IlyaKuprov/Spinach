% MEX compilation utility. Rebuilds all C++ and CUDA MEX binaries of
% the Spinach installation that contains this file, placing
% each binary next to its source. CUDA compilation requires
% mexcuda and a supported CUDA toolkit; on Windows, CUDA_PATH
% must point to that toolkit. Syntax:
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
mex('-R2018a','-O','-DNDEBUG','COMPFLAGS=$COMPFLAGS','LINKFLAGS=$LINKFLAGS',...
    [P '/kernel/line_shapes/lorentzcon.cpp'],'-outdir',[P '/kernel/line_shapes']);

% Gaussian convolution
mex('-R2018a','-O','-DNDEBUG','COMPFLAGS=$COMPFLAGS','LINKFLAGS=$LINKFLAGS',...
    [P '/kernel/line_shapes/gausscon.cpp'],'-outdir',[P '/kernel/line_shapes']);

% Cubic polynomial roots
mex('-R2018a','-O','-DNDEBUG','COMPFLAGS=$COMPFLAGS','LINKFLAGS=$LINKFLAGS',...
    [P '/kernel/eigenfields/cubic_roots.cpp'],'-outdir',[P '/kernel/eigenfields']);

% Sparse double row sorter
mex('-R2018a','-O','-DNDEBUG',...
    [P '/kernel/indexing/spsortrows.cpp'],'-outdir',[P '/kernel/indexing']);

% Sparse double unique columns
mex('-R2018a','-O','-DNDEBUG',...
    [P '/kernel/indexing/spunicols.cpp'],'-outdir',[P '/kernel/indexing']);

% CUDA sparse matrix product
if ispc
    mexcuda('-R2018a','-O','-lcuda',...
            ['-L' fullfile(getenv('CUDA_PATH'),'lib','x64')],...
            [P '/kernel/overloads/cuda_sparse_by_sparse_mex.cu'],...
            '-outdir',[P '/kernel/overloads']);
else
    mexcuda('-R2018a','-O','-lcuda',...
            [P '/kernel/overloads/cuda_sparse_by_sparse_mex.cu'],...
            '-outdir',[P '/kernel/overloads']);
end

end

% Audiophiles don't use their equip-
% ment to listen to your music. Audio-
% philes use your music to listen to
% their equipment.
%
% Alan Parsons


