% This function runs contiouously on one of our servers
% at Southampton, and its purpose is to catch any unin-
% tended consequences before they propagate too far down
% the development chain. Syntax:
%
%                  patrol(test_subject)
%
% Parameters:
%
%    test_subject - a character string; if it occurs
%                   anywhere within the example file
%                   path, that file is included into
%                   the patrol run
%
% Outputs:
%
%    whatever the individual examples return
% 
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=patrol.m>

function patrol(test_subject)

% Set default and check consistency
if ~exist('test_subject','var')
    test_subject='';
end
grumble(test_subject);

% List exceptions
exceptions={};

% Shuffle the RNG
rng('shuffle');

% Get the directory tree
P=mfilename('fullpath'); P=P(1:(end-6));
mfiles=dir([P '../../examples/**/*.m']);

% Find relevant files
if ~isempty(test_subject)
    
    % No file is initially relevant
    relevant_file_mask=false(numel(mfiles),1);
    
    % Loop over the file list
    for n=1:numel(mfiles)
        
        % Get full name
        full_name=[mfiles(n).folder filesep mfiles(n).name];
        
        % Read the file
        fid=fopen(full_name,'r');
        content=textscan(fid,'%s','Delimiter','\n');
        fclose(fid); content=content{1};
        
        % Files that mention the subject are relevant
        for m=1:numel(content)
            if contains(content{m},test_subject)||...
               contains(full_name,test_subject)
                relevant_file_mask(n)=true();
            end
        end
        
    end
    
else
    
    % Test all example files
    relevant_file_mask=true(numel(mfiles),1);
    
end

% Update the file list
mfiles=mfiles(relevant_file_mask);

% Inform the user
if isempty(mfiles)
    error('the file list is empty.');
else
    disp(['Patrol over ' num2str(numel(mfiles)) ' files...']);
end

% Enter the main loop
hashes_match=true();
while hashes_match
    
    % Pick a random file
    n=randi(numel(mfiles));

    % Check that Matlab's syntax checker is on green
    if ~isempty(checkcode([mfiles(n).folder filesep mfiles(n).name]))
        edit(file_name); error('the built-in syntax checker has something to say');
    end
    
    % Run the file
    if ~ismember(mfiles(n).name,exceptions)
        disp(['Running ' mfiles(n).name ' ...']);
        cd(mfiles(n).folder); 
        eval(mfiles(n).name(1:(end-2)));
    end
    
    % Flush the display buffer
    drawnow(); pause(1);
    
end

end

% Consistency enforcement
function grumble(test_subject)
if ~ischar(test_subject)
    error('test_subject must be a character string');
end
end

% No snowflake in an avalanche ever feels responsible.
%
% Stanislav Lec

