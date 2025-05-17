% Copies all m-files from everywhere in Spinach into a 
% specified directory; this is useful for getting LLMs
% to review the entire code base. Syntax:
%
%                   gpt_pack(dir_name)
%
% Parameters:
%
%    dir_name - directory name
%
% Ouputs:
%
%    files copied to the specified directory
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=gpt_pack.m>

function gpt_pack(dir_name)

% Check consistency
grumble(dir_name);

% Get the directory tree
P=mfilename('fullpath'); P=P(1:(end-9));
mfiles=dir([P '..' filesep '..' filesep ...
              '..' filesep '**' filesep '*.m']);

% Relative path start
rps_cut=numel(P)-17;

% Create directories
for n=1:numel(mfiles)
    folder=[dir_name filesep mfiles(n).folder(rps_cut:end)];
    if ~exist(folder,'dir'); mkdir(folder); end
end

% Copy files
for n=1:numel(mfiles)
     source_file=[mfiles(n).folder filesep mfiles(n).name];
     destin_file=[dir_name filesep mfiles(n).folder(rps_cut:end) ...
                  filesep mfiles(n).name];
     copyfile(source_file,destin_file);
end

end

% Consistency enforcement
function grumble(dir_name)
if ~ischar(dir_name)
    error('dir_name must be a character string.');
end
end

% Sheep spend their lives being afraid of wolves,
% and get eaten by shepherds.
%
% Georgian saying

