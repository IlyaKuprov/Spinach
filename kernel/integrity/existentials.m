% Kernel integrity control. Checks for collisions between Spinach 
% functions and anything else that the user may have installed or
% written in the current Matlab instance. Also checks for any fi-
% les that are not visible to Matlab because the corresponding di-
% rectory is not on the path.
%
% Collisions of function names and path problems are the most fre-
% quent support topic at the forum.
%
% ilya.kuprov@weizmann.ac.il
% david.goodwin@kit.edu
%
% <https://spindynamics.org/wiki/index.php?title=existentials.m>

function existentials()

% Do not run inside parallel pools
if isworkernode, return; end

% Inform the user
disp('Running startup checks...');

%##########################################
% NO, IT WILL NOT MAGICALLY START WORKING %
%     IF YOU COMMENT ANY OF THIS OUT      %
%##########################################

% Global Matlab release
if isMATLABReleaseOlderThan('R2024b')
    error('Spinach requires Matlab R2024b or later.');
end

% Existential toolboxes
if ~exist([matlabroot filesep 'toolbox' filesep 'parallel'],'dir')
    error('Spinach requires Parallel Computing Toolbox.');
end
if ~exist([matlabroot filesep 'toolbox' filesep 'nnet'],'dir')
    error('Spinach requires Deep Learning Toolbox.');
end
if ~exist([matlabroot filesep 'toolbox' filesep 'rl'],'dir')
    error('Spinach requires Reinforcement Learning Toolbox.');
end
if ~exist([matlabroot filesep 'toolbox' filesep 'optim'],'dir')
    error('Spinach requires Optimisation Toolbox.');
end
if ~exist([matlabroot filesep 'toolbox' filesep 'stats'],'dir')
    error('Spinach requires Statistics and Machine Learning Toolbox.');
end
if ~exist([matlabroot filesep 'toolbox' filesep 'map'],'dir')
    error('Spinach requires Mapping Toolbox.');
end
if ~exist([matlabroot filesep 'toolbox' filesep 'aero'],'dir')
    error('Spinach requires Aerospace Toolbox.');
end

% List top level directories
top_level={'kernel','interfaces','experiments','etc'};

% Get the directory trees
mfiles=cell(numel(top_level,1));
P=mfilename('fullpath'); P=P(1:(end-13));
for n=1:numel(top_level)
    mfiles{n}=dir([P filesep '..' filesep '..' filesep ...
                   top_level{n} filesep '**' filesep '*.m']);
end

% Process the files
for n=1:numel(top_level)
    for k=1:numel(mfiles{n})
        
        % Get the full name as per Spinach
        file_name_spinach=[mfiles{n}(k).folder filesep mfiles{n}(k).name];
        
        % Get the full name as per Matlab
        file_name_matlab=which(mfiles{n}(k).name);
        
        % Check for collisions, except in overloads
        if (~isempty(file_name_matlab))&&...
           (~strcmp(file_name_spinach,file_name_matlab))&&...
           (~contains(mfiles{n}(k).folder,'overloads'))
            disp('Function name collision problem.');
            disp('============');
            disp(['File: ' mfiles{n}(k).name]);
            disp(['Expected place: ' file_name_spinach]);
            disp(['Matlab reports: ' file_name_matlab]);
            disp('============');
            disp('Please remove or rename the external file.');
            error('startup checks not passed');
        end
        
        % Check for missing files
        if isempty(file_name_matlab)
            disp('Matlab path setup problem.');
            disp('============');
            disp(['File: ' mfiles{n}(k).name]);
            disp(['Expected place: ' file_name_spinach]);
            disp('Matlab says the file is not in its path.');
            disp('============');
            disp('Please make sure the path is correctly set.');
            error('startup checks not passed');
        end
        
    end
end

% Windows
if ispc

    % Find out which file system Spinach volume has
    own_disk=mfilename('fullpath'); own_disk=own_disk(1:3);
    own_disk=System.IO.DriveInfo(own_disk);

    % Warn about non-NTFS volumes
    if ~strcmp(char(own_disk.DriveFormat),'NTFS')
        warning('Spinach:NonNTFSVolume',['Spinach is running from a ' ...
                char(own_disk.DriveFormat) ' volume; NTFS is recommended ' ...
                'on Windows for reliable file locking and performance.']);
    end

end

end

% It has been my observation that most people get ahead
% during the time that others waste.
%
% Henry Ford

% #NHEAD #NGRUM

