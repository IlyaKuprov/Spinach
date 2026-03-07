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
    mfiles{n}=dir([P '..' filesep '..' filesep ...
                   top_level{n} filesep '**' filesep '*.m']);
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
        raw_content=textscan(fid,'%s','Delimiter','\n','Whitespace','');
        fclose(fid); raw_content=raw_content{1};
        content=deblank(raw_content);

        % Identify top-level construct
        code_line_idx=find((~cellfun(@isempty,content))&...
                           (~startsWith(strtrim(content),'%')),1,'first');
        if isempty(code_line_idx)
            top_line='';
            is_function=false;
            is_classdef=false;
            has_input_args=false;
            has_output_args=false;
        else
            top_line=strtrim(content{code_line_idx});
            is_function=startsWith(top_line,'function');
            is_classdef=startsWith(top_line,'classdef');
            arg_block=regexp(top_line,'\(([^)]*)\)','tokens','once');
            has_input_args=(~isempty(arg_block))&&(~isempty(strtrim(arg_block{1})));
            has_output_args=contains(top_line,'=');
        end

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

        % Check for hard tab characters
        for m=1:numel(raw_content)
            if contains(raw_content{m},sprintf('\t'))
                edit(file_name); error(['tab character found in line ' num2str(m)]);
            end
        end

        % Check for two line breaks at end of file
        file_text=fileread(file_name);
        if isempty(regexp(file_text,'(\r\n|\n){2}$','once'))
            edit(file_name); error('file must end with two blank lines');
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

        % Check that grumbler call is present in main body
        if (~grumbler_ex)&&is_function&&(~is_classdef)&&has_input_args&&...
           (~contains(mfiles(k).folder,{'legacy','includes','overloads'}))

            % Locate local function declarations
            function_lines=find(startsWith(strtrim(content),'function'));

            % Find end of main function
            if numel(function_lines)>1
                main_fun_end=function_lines(2)-1;
            else
                main_fun_end=numel(content);
            end

            % Find grumble call in main function body
            grumble_call=false;
            for m=(code_line_idx+1):main_fun_end
                if startsWith(strtrim(content{m}),'grumble(')
                    grumble_call=true;
                    break;
                end
            end
            if ~grumble_call
                edit(file_name); error('grumbler call missing from main function body');
            end

            % Ensure grumble is the last function in file
            grumble_line=find(contains(content,'function grumble('),1,'first');
            if ~isempty(grumble_line)
                if function_lines(end)~=grumble_line
                    edit(file_name); error('grumbler function must be last in file');
                end

                % Find end of grumble function
                grumble_end=find_block_end(content,grumble_line);
                if isempty(grumble_end)
                    edit(file_name); error('unable to determine grumbler function end');
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
            pause(0.25); % do not make too many requests all at once
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
        if ~header_ex

            % Locate end of documentation header
            doc_end=find(cellfun(@isempty,content),1,'first');
            if isempty(doc_end)
                edit(file_name); error('documentation header must end with a blank line');
            end
            doc_header=content(1:(doc_end-1));

            % Require syntax section
            if ~any(contains(doc_header,'Syntax:'))
                edit(file_name); error('documentation header must contain Syntax section');
            end

            % Require parameters section when inputs are present
            if is_function&&has_input_args&&...
              (~any(contains(doc_header,'Parameters:')))&&...
              (~any(contains(doc_header,'Arguments:')))
                edit(file_name); error('documentation header must contain Parameters or Arguments section');
            end

            % Require output section when outputs are present
            if is_function&&has_output_args&&...
              (~any(contains(doc_header,'Output:')))&&...
              (~any(contains(doc_header,'Outputs:')))
                edit(file_name); error('documentation header must contain Output or Outputs section');
            end

        end
        
        % Enforce filesep in all paths
        path_calls={'dir(','fopen(','load(','save(','copyfile(',...
                    'movefile(','delete(','mkdir(','rmdir(',...
                    'cd(','addpath(','genpath(','exist('};
        for m=1:numel(content)

            % Skip blank lines and comments
            line_text=strtrim(content{m});
            if isempty(line_text), continue; end
            if strcmp(line_text(1),'%'), continue; end

            % Process file operation call sites only
            if ~any(contains(line_text,path_calls)), continue; end

            % Process string literals in the line
            string_literals=regexp(line_text,'''[^'']*''','match');
            for n=1:numel(string_literals)
                token=string_literals{n}(2:(end-1));
                if contains(token,'http://')||contains(token,'https://')
                    continue;
                end
                if contains(token,'/')||contains(token,'\')
                    edit(file_name); error(['hard-coded path separator in line ' num2str(m)]);
                end
            end

        end
        
        % Otherwise clause in each switch block
        for m=1:numel(content)

            % Check each switch declaration
            line_text=strtrim(content{m});
            if isempty(regexp(line_text,'^switch\b','once'))
                continue;
            end

            % Parse the switch block
            block_depth=1;
            has_otherwise=false;
            for n=(m+1):numel(content)

                % Tokenise current line
                [tokens,token_start,token_end,code_line]=control_tokens(content{n});

                % Process tokens in lexical order
                for p=1:numel(tokens)

                    % Block-opening keywords
                    if is_block_start_token(tokens{p})
                        block_depth=block_depth+1;
                        continue;
                    end

                    % Top-level otherwise clause
                    if strcmp(tokens{p},'otherwise')&&(block_depth==1)
                        has_otherwise=true;
                        continue;
                    end

                    % Block-closing end keyword
                    if strcmp(tokens{p},'end')&&...
                       is_block_end_token(code_line,token_start(p),token_end(p))
                        block_depth=block_depth-1;
                        if block_depth==0
                            break;
                        end
                    end

                end

                % Exit when switch block terminates
                if block_depth==0
                    break;
                end

            end

            % Reject malformed switch blocks
            if block_depth~=0
                edit(file_name); error(['unterminated switch block at line ' num2str(m)]);
            end

            % Require otherwise at top level
            if ~has_otherwise
                edit(file_name); error(['switch block without otherwise at line ' num2str(m)]);
            end

        end
        
        % disp() instead of report() when spin_system is available
        if is_function&&contains(top_line,'spin_system')
            for m=1:numel(content)
                line_text=strtrim(content{m});
                if isempty(line_text), continue; end
                if strcmp(line_text(1),'%'), continue; end
                if contains(line_text,'disp(')
                    edit(file_name); error(['disp() used when spin_system is available, line ' num2str(m)]);
                end
            end
        end
        
    end
end

% Report success
disp('Now ye are clean through the word which I have spoken unto you. - John 15:3');

end

% Detect block-opening MATLAB keywords
function answer=is_block_start_token(token)
answer=ismember(token,{'if','for','parfor','while',...
                       'switch','try','spmd','function','classdef'});
end

% Extract control-flow tokens from a line
function [tokens,token_start,token_end,code_line]=control_tokens(line_text)
code_line=lower(strip_strings_and_comments(line_text));
[token_start,token_end,tokens]=regexp(code_line,...
    '(?<![A-Za-z0-9_])(if|for|parfor|while|switch|try|spmd|function|classdef|otherwise|end)(?![A-Za-z0-9_])',...
    'start','end','match');
end

% Remove string literals and comments from a line
function code_line=strip_strings_and_comments(line_text)
code_chars=[]; in_string=false; n=1;
while n<=numel(line_text)
    if in_string
        if strcmp(line_text(n),'''')
            if (n<numel(line_text))&&strcmp(line_text(n+1),'''')
                n=n+2;
            else
                in_string=false;
                n=n+1;
            end
        else
            n=n+1;
        end
    else
        if strcmp(line_text(n),'''')
            in_string=true;
            n=n+1;
        elseif strcmp(line_text(n),'%')
            break;
        else
            code_chars=[code_chars line_text(n)]; %#ok<AGROW>
            n=n+1;
        end
    end
end
code_line=char(code_chars);
end

% Decide whether an "end" token closes a code block
function answer=is_block_end_token(code_line,start_idx,end_idx)
prev_idx=find(~isspace(code_line(1:(start_idx-1))),1,'last');
next_rel=find(~isspace(code_line((end_idx+1):end)),1,'first');
if isempty(prev_idx)
    prev_char=[];
else
    prev_char=code_line(prev_idx);
end
if isempty(next_rel)
    next_char=[];
else
    next_char=code_line(end_idx+next_rel);
end
answer=(isempty(prev_char)||any(prev_char==[';',',']))&&...
       (isempty(next_char)||any(next_char==[';',',']));
end

% Find terminating end for a block
function end_line=find_block_end(content,start_line)
block_depth=1; end_line=[];
for n=(start_line+1):numel(content)

    % Tokenise current line
    [tokens,token_start,token_end,code_line]=control_tokens(content{n});

    % Process control-flow tokens
    for p=1:numel(tokens)
        if is_block_start_token(tokens{p})
            block_depth=block_depth+1;
        elseif strcmp(tokens{p},'end')&&...
               is_block_end_token(code_line,token_start(p),token_end(p))
            block_depth=block_depth-1;
            if block_depth==0
                end_line=n;
                return;
            end
        end
    end

end
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
