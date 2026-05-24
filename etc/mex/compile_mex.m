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
here=fileparts(mfilename('fullpath'));
spinach_root=fileparts(fileparts(here));

% Get source and output locations
src_files={fullfile(spinach_root,'kernel','line_shapes','lorentzcon.cpp'),...
           fullfile(spinach_root,'kernel','eigenfields','cubic_roots.cpp')};
out_dirs={fullfile(spinach_root,'kernel','line_shapes'),...
          fullfile(spinach_root,'kernel','eigenfields')};

% Check source availability
src_missing=~cellfun(@isfile,src_files);
if all(src_missing)
    error('No C++ source files found.');
elseif any(src_missing)
    error('Required C++ source file is missing.');
end

% Compile with case-specific options
mex_args={'-R2018a','-O','COMPFLAGS=$COMPFLAGS','LINKFLAGS=$LINKFLAGS'};
for n=1:numel(src_files)
    mex(mex_args{:},src_files{n},'-outdir',out_dirs{n});
end

end

% Audiophiles don't use their equip-
% ment to listen to your music. Audio-
% philes use your music to listen to
% their equipment.
%
% Alan Parsons
