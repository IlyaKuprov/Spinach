% MEX compilation utility. Rebuilds all C++ MEX binaries in the current
% directory. Syntax:
%
%                      compile_mex()
%
% Parameters:
%
%   none
%
% Outputs:
%
%   none

function compile_mex()

% Get C++ source file list
cpp_files=dir('*.cpp');

% Halt if there is nothing to compile
if isempty(cpp_files)
    error('No C++ source files were found in the current directory.');
end

% Rebuild all MEX binaries in the current directory
for n=1:numel(cpp_files)
    mex('-R2018a','-O',cpp_files(n).name,'-outdir',pwd);
end

end


