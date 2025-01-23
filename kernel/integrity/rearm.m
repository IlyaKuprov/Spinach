% Rearms the sniffer database. The sniffer checks Spinach 
% distribution .m files for any modifications that the user
% did since downloading Spinach. The function prints the 
% list of files that have changed in any way since the in-
% ternal database has been rearmed.
%
% The purpose is to catch local modifications that the user
% may have made and forgotten about, that are causing some un-
% intended consequences elsewhere in Spinach.
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=rearm.m>

function rearm()

% List top level directories
top_level={'kernel','interfaces','experiments','etc'};

% List exceptions
exceptions={};

% Get the directory trees
mfiles=cell(numel(top_level,1));
P=mfilename('fullpath'); P=P(1:(end-5));
for n=1:numel(top_level)
    mfiles{n}=dir([P '../../' top_level{n} '/**/*.m']);
end

% Get the table going
smells={};

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
            smells{end+1}=md5_hash([mfiles{n}(k).name md5_hash(content)]); %#ok<AGROW>
            
        end
    end
end

% Overwrite the file
delete([P 'smells.mat']); save([P 'smells.mat'],'smells');

% Inform the user
disp('rearm: sniffer rearmed.');

end

% MEN WANTED for hazardous journey, small wages, bitter cold, long months
% of complete darkness, constant danger, safe return doubtful, honour and
% recognition in case of success.
%
% Ernest Shackleton's advertisement in The Times

