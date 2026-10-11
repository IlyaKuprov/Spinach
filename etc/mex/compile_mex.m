% MEX compilation utility. Rebuilds all C++ and CUDA MEX binaries of
% the Spinach installation that contains this file, placing
% each binary next to its source. CUDA compilation requires
% MATLAB R2026b or later, mexcuda, a supported host C++ compiler,
% and the CUDA driver library. Windows uses CUDA_PATH/lib/x64 when
% CUDA_PATH is set; otherwise mexcuda's library search is used.
% macOS CUDA-unavailable errors are reported after the C++ builds;
% all other compilation errors propagate.
% Syntax:
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

% Enforce the CUDA sparse-storage release contract
if isMATLABReleaseOlderThan('R2026b')
    disp('CUDA sparse-product MEX compilation skipped: requires R2026b or later.');
    return;
end

% CUDA driver library search
lib_args={};
if ispc&&~isempty(getenv('CUDA_PATH'))
    lib_args={['-L' fullfile(getenv('CUDA_PATH'),'lib','x64')]};
end

% CUDA sparse matrix product
try
    mexcuda('-R2018a','-O',lib_args{:},'-lcuda',...
            [P '/kernel/overloads/cuda_sparse_by_sparse_mex.cu'],...
            '-outdir',[P '/kernel/overloads']);
catch exception
    if ismac&&strcmp(exception.identifier,'parallel:gpu:device:MacOSGPUNotSupported')
        disp(['CUDA sparse-product MEX unavailable: ' exception.message]);
    else
        rethrow(exception);
    end
end

end

% Audiophiles don't use their equip-
% ment to listen to your music. Audio-
% philes use your music to listen to
% their equipment.
%
% Alan Parsons


