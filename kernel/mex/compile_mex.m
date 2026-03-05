% MEX compilation utility. Rebuilds all C++ MEX binaries in 
% the current directory. Syntax:
%
%                       compile_mex()
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=compile_mex.m>

function compile_mex() % #NGRUM #NHEAD

% Build compiler argument list
if ispc
    mex_args={'-R2018a','-O',...
              'COMPFLAGS=$COMPFLAGS /openmp',...
              'LINKFLAGS=$LINKFLAGS /openmp'};
elseif isunix&&(~ismac)
    mex_args={'-R2018a','-O',...
              'CXXFLAGS=$CXXFLAGS -fopenmp',...
              'LDFLAGS=$LDFLAGS -fopenmp'};
else
    mex_args={'-R2018a','-O'};
end

% Get this function directory
here=fileparts(mfilename('fullpath'));

% Get C++ source file list in this function directory
cpp_files=dir(fullfile(here,'*.cpp'));

% Halt if there is nothing to compile
if isempty(cpp_files)
    error('No C++ source files were found in %s.',here);
end

% Rebuild all MEX binaries in this function directory
for n=1:numel(cpp_files)
    src_file=fullfile(here,cpp_files(n).name);
    mex(mex_args{:},src_file,'-outdir',here);
end

end

% Audiophiles don't use their equip-
% ment to listen to your music. Audio-
% philes use your music to listen to
% their equipment.
%
% Alan Parsons

