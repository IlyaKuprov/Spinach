% Kernel integrity control. Checks Spinach distribution .m files
% for any modifications that the user did since downloading Spi-
% nach. The function prints the list of files that have changed
% in any way since the internal database has been rearmed.
%
% The purpose is to catch local modifications that the user may
% have made and forgotten about, that are causing some unintend-
% ed consequences elsewhere in Spinach.
%
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=sniff.m>

function sniff(action)

% Default is to take no action
if ~exist('action','var'), action='none'; end

% Check consistency
grumble(action);

% List top level directories
top_level={'kernel','interfaces','experiments','etc'};

% List exceptions
exceptions={};

% Get the directory trees
mfiles=cell(numel(top_level,1));
P=mfilename('fullpath'); P=P(1:(end-5));
for n=1:numel(top_level)
    mfiles{n}=dir([P filesep '..' filesep '..' filesep ...
                   top_level{n} filesep '**' filesep '*.m']);
end

% Read the smells table
load('smells.mat','smells'); all_good=true();

% Start line counters
code_line_count=0;
comm_line_count=0;

% Process the files
for n=1:numel(top_level)
    for k=1:numel(mfiles{n})
        if ~ismember(mfiles{n}(k).name,exceptions)
            
            % Read the file
            file_name=[mfiles{n}(k).folder filesep mfiles{n}(k).name];
            fid=fopen(file_name,'r');
            content=textscan(fid,'%s','Delimiter', '\n');
            fclose(fid); content=content{1};
            
            % Get the smell
            smell=md5_hash([mfiles{n}(k).name md5_hash(content)]);

            % Drop blank lines
            content(cellfun(@(x)isempty(deblank(x)),content))=[];
            
            % Count comment lines
            comm_line_count=comm_line_count+...
                            nnz(cellfun(@(x)strcmp('%',x(1)),content));

            % Drop comment lines
            content(cellfun(@(x)strcmp('%',x(1)),content))=[];

            % Count code lines
            code_line_count=code_line_count+numel(content);
            
            % Check for a match
            if ~ismember(smell,smells)
                
                % Previously unseen file
                disp(['smells fishy: ' file_name]); all_good=false();
                
                % Open the file
                if strcmp(action,'open'), edit(file_name); end
                    
            end
            
        end
    end
end

% Give an all-clear if appropriate
if all_good
    disp(['sniff: ' num2str(comm_line_count) ' lines of comments,']);
    disp(['sniff: ' num2str(code_line_count) ' lines of code,']);
    disp('sniff: everything smells fine.');
end

end

% Consistency enforcement
function grumble(action)
if (~ischar(action))||(~ismember(action,{'none','open'}))
    error('action must be ''none'' or ''open''.');
end
end

% Glen had an unusual introduction to the Arctic in 1932. He
% thought he had accepted a friend's invitation to a debutante
% dance, then found that it was to go to Spitzbergen as one of
% the eight-man crew of a 45-foot Peterhead fishing boat owned
% by a Cambridge law don. The expedition committed him to 4000
% miles of sailing and two months of surveying in the mounta-
% ins; it left him fascinated by the Arctic.
%
% The obituary to Sir Alexander Glen,
% published by The Telegraph.

