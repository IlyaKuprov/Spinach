function compile_prune_cpu()
%COMPILE_PRUNE_CPU Build prune_cpu MEX for the current platform.
%
%   compile_prune_cpu

here = fileparts(mfilename('fullpath'));
src  = fullfile(here,'prune_cpu.cpp');

if ~exist(src,'file')
    error('Source file not found: %s', src);
end

% Use the R2017b MEX API to allow mxUnshareArray() (copy-on-write detach).
% This targets the classic (separate real/imag) complex API, which is fine for Spinach.
mex('-R2017b','-O',src,'-outdir',here);

end
