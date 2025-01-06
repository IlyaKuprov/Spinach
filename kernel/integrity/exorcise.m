% Searches Spinach distribution folders for any functions that
% do not conform to the house style. Opens the first one and
% complains to the console. Syntax:
%
%                        exorcise(mode)
%
% Parameters:
%
%   mode   - 'online' (default) checks the documentation
%            Wiki for the corresponding page; 'offline'
%            skips the Wiki check 
%
% Any user contribution that this function has something to
% say about will either be brought under the house style, or
% rejected back to the user, depending on the amount of work
% involved. Always run this function before a commit if you
% have write access to Spinach repository.
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=exorcise.m>

function exorcise(mode)

% Default mode is online
if ~exist('mode','var'), mode='online'; end

% Check consistency
grumble(mode);

% List top level directories
top_level={'kernel','interfaces','experiments','etc'};

% External package directories
foreign_dirs={'jsonlab-1.5'};

% Get the directory trees
mfiles=cell(numel(top_level,1));
P=mfilename('fullpath'); P=P(1:(end-8));
for n=1:numel(top_level)
    mfiles{n}=dir([P '../../' top_level{n} '/**/*.m']);
end
mfiles=cell2mat(mfiles'); 
mfiles=mfiles(randperm(numel(mfiles)));

% Process the files
for k=1:numel(mfiles)
    
    % Update the user
    disp(['inspecting file ' num2str(k) ' out of ' num2str(numel(mfiles)) '...']);
    
    % Exclude foreign package directories
    if ~any(cellfun(@(x)contains(mfiles(k).folder,x),foreign_dirs))
        
        % Read the file
        file_name=[mfiles(k).folder filesep mfiles(k).name];
        fid=fopen(file_name,'r');
        content=textscan(fid,'%s','Delimiter', '\n');
        fclose(fid); content=deblank(content{1});

        % Apply exceptions
        grumbler_ex=false; header_ex=false; wiki_ex=false;
        for m=1:numel(content)
            if contains(content{m},'#NGRUM')
                grumbler_ex=true;
            end
            if contains(content{m},'#NHEAD')
                header_ex=true;
            end
            if contains(content{m},'#NWIKI')
                wiki_ex=true;
            end
        end
        
        % Check for excessive white space
        blank_lines=cellfun(@isempty,content);
        for m=2:numel(blank_lines)
            if blank_lines(m-1)&&blank_lines(m)
                edit(file_name); error(['two sequential blank lines at ' num2str(m)]);
            end
        end
        
        % Check that the grumbler is present
        if ~grumbler_ex
            for m=1:numel(content)
                if contains(content{m},'function grumble(')||...
                   contains(content{m},'classdef')||...
                   contains(content{m},'% Catch incorrect calls')||...
                  (contains(content{m},'function')&&contains(content{m},'()'))||...
                   contains(mfiles(k).folder,{'legacy','includes','overloads'})
                    break;
                end
                if m==numel(content)
                    edit(file_name); error('grumbler function missing');
                end
            end
        end
        
        % Enforce the minimum doc header size
        if ~header_ex
            doc_header_lines=0;
            for m=1:numel(content)
                if (doc_header_lines<=9)&&(isempty(content{m})||(~strcmp(content{m}(1),'%')))
                    edit(file_name); error('documentation header too short');
                else
                    doc_header_lines=doc_header_lines+1;
                end
                if doc_header_lines>9, break; end
            end
        end
        
        % Check that function header has a Wiki link
        if ~wiki_ex
            for m=1:numel(content)
                if isempty(content{m})
                    if (numel(content{m-1})<7)||(~strcmp(content{m-1}(1:7),'% <http'))
                        edit(file_name); error('last line of the header must be the wiki link');
                    else
                        wiki_address=content{m-1}(4:(end-1)); break;
                    end
                end
            end
        end
        
        % Check that Matlab's syntax checker is on green
        if ~isempty(checkcode(file_name))
            edit(file_name); error('the built-in syntax checker has something to say');
        end
        
        % Check that the documentation page exists on the Wiki
        if strcmp(mode,'online')&&(~wiki_ex)
            connection=java.net.URL(wiki_address).openConnection();
            connection.setDoInput(true); connection.setDoOutput(true);
            connection.setRequestProperty('User-Agent','spinach_exorcism_module');
            try
                input=java.io.BufferedReader(...
                      java.io.InputStreamReader(connection.getInputStream()));
                input.close();
            catch
                edit(file_name); web(wiki_address,'-browser');
                error('documentation page missing from the wiki');
            end
            pause(0.5); % do not make too many requests all at once
        end
        
        % Check that norms specify type explicitly
        for m=1:numel(content)
            if (~isempty(content{m}))&&contains(content{m},'norm(')&&...               % #NORMOK
                 isempty(regexp(content{m},'[norm(]*[,]\d[)]','once'))&&...            % #NORMOK
               (~contains(content{m},'fro'))&&(~contains(content{m},'inf'))&&...
               (~contains(content{m},'cheap'))&&(~contains(content{m},'#NORMOK'))
                edit(file_name); error(['unspecified norm type in line ' num2str(m)]);
            end
        end
           
        % Check keywords in the docs header
        
        % Spaces around arithmetical and logical operations
        
        % Brackets around intervals
        
        % Overlong lines
        
        % Enforce filesep in all paths
        
        % Otherwise clause in each switch block
        
        % Quotation present
        
        % Percentage sign with no space afterwards, or multiple percentage signs
        
        % disp() instead of report() when spin_system is available
        
        % RNG calls
        
        % No white line before a comment
        
    end
end

% Report success
disp('Now ye are clean through the word which I have spoken unto you. - John 15:3');

end

% Consistency enforcement
function grumble(mode)
if ~ischar(mode)
    error('mode must be a character string.');
end
if ~ismember(mode,{'online','offline'})
    error('mode can be either ''online'' or ''offline''.');
end
end

% It is my ambition to say in ten sentences what everyone 
% else says in a whole book - what everyone else does not
% say in a whole book.
%
% Friedrich Nietzsche

