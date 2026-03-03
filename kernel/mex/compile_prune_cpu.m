function compile_prune_cpu()
%COMPILE_PRUNE_CPU Build prune_cpu MEX for the current platform.
%
%   compile_prune_cpu

here = fileparts(mfilename('fullpath'));
src  = fullfile(here,'prune_cpu.cpp');

if ~exist(src,'file')
    error('Source file not found: %s', src);
end

% Use modern interleaved-complex MEX API (R2018a+).
% Spinach minimum supported MATLAB is R2024b.
mex('-R2018a','-O',src,'-outdir',here);

end
