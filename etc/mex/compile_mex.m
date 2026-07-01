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

% Compile sparse row sorter with platform OpenMP flags
if ispc

    % Query the selected Windows C++ compiler
    cc=mex.getCompilerConfigurations('C++','Selected');
    cc_name=lower([cc.Name ' ' cc.Manufacturer]);

    % Use GCC-compatible flags for MinGW
    if contains(cc_name,'mingw')||contains(cc_name,'gcc')
        mex('-R2018a','-O','-DNDEBUG','CXXFLAGS=$CXXFLAGS -fopenmp','LDFLAGS=$LDFLAGS -fopenmp',...
            [P '/kernel/indexing/spsortrows.cpp'],'-outdir',[P '/kernel/indexing']);

    % Use MSVC OpenMP flags otherwise
    else
        mex('-R2018a','-O','-DNDEBUG','COMPFLAGS=$COMPFLAGS /openmp',...
            [P '/kernel/indexing/spsortrows.cpp'],'-outdir',[P '/kernel/indexing']);
    end

elseif ismac

    % Locate Homebrew or MacPorts OpenMP runtime
    omp_inc=''; omp_lib='';
    if isfile('/opt/homebrew/opt/libomp/include/omp.h')
        omp_inc='/opt/homebrew/opt/libomp/include';
        omp_lib='/opt/homebrew/opt/libomp/lib';
    elseif isfile('/usr/local/opt/libomp/include/omp.h')
        omp_inc='/usr/local/opt/libomp/include';
        omp_lib='/usr/local/opt/libomp/lib';
    elseif isfile('/opt/homebrew/include/omp.h')
        omp_inc='/opt/homebrew/include';
        omp_lib='/opt/homebrew/lib';
    elseif isfile('/usr/local/include/omp.h')
        omp_inc='/usr/local/include';
        omp_lib='/usr/local/lib';
    elseif isfile('/opt/local/include/libomp/omp.h')
        omp_inc='/opt/local/include/libomp';
        omp_lib='/opt/local/lib/libomp';
    elseif isfile('/opt/local/include/omp.h')
        omp_inc='/opt/local/include';
        omp_lib='/opt/local/lib';
    end

    % Fall back to serial compilation if OpenMP is missing
    if isempty(omp_inc)
        warning('compile_mex:openmpMissing',...
                'OpenMP libomp was not found on macOS; compiling spsortrows without OpenMP.');
        mex('-R2018a','-O','-DNDEBUG',...
            [P '/kernel/indexing/spsortrows.cpp'],'-outdir',[P '/kernel/indexing']);

    % Use Apple Clang OpenMP flags otherwise
    else
        mex('-R2018a','-O','-DNDEBUG',...
            ['CXXFLAGS=$CXXFLAGS -Xpreprocessor -fopenmp -I' omp_inc],...
            ['LDFLAGS=$LDFLAGS -L' omp_lib ' -Wl,-rpath,' omp_lib ' -lomp'],...
            [P '/kernel/indexing/spsortrows.cpp'],'-outdir',[P '/kernel/indexing']);
    end

elseif isunix

    % Use GCC or Clang OpenMP flags
    mex('-R2018a','-O','-DNDEBUG','CXXFLAGS=$CXXFLAGS -fopenmp','LDFLAGS=$LDFLAGS -fopenmp',...
        [P '/kernel/indexing/spsortrows.cpp'],'-outdir',[P '/kernel/indexing']);

else

    % Fall back to serial compilation on unknown platforms
    mex('-R2018a','-O','-DNDEBUG',...
        [P '/kernel/indexing/spsortrows.cpp'],'-outdir',[P '/kernel/indexing']);
end

end

% Audiophiles don't use their equip-
% ment to listen to your music. Audio-
% philes use your music to listen to
% their equipment.
%
% Alan Parsons
