% Absence of the retired CCE interface and its executable callers.
% Syntax: result=test_cce_no_legacy()
%
% Checks the active path and Git-tracked MATLAB sources below the repository
% containing this test. Comments and quoted text are stripped before the
% direct-call and function-handle scan; transpose operators remain code.
%
% ilya.kuprov@weizmann.ac.il

function result=test_cce_no_legacy()

% Check that no removed entry point is visible on the active path
result=new_test_result('kernel/cce_no_legacy','Native CCE interface only',...
                       'Retired names, tracked callers, and native entry points.');
retired={'cce','cce_bath','cce_cluster','cce_master','cce_noise',...
         'cce_observe','cce_prepare','cce_sequence','cce_validate'};
for n=1:numel(retired)
    result=test_true(result,['Removed ' retired{n}],...
                     ismember(exist(retired{n},'file'),[0 7])&&isempty(which(retired{n})),...
                     'The retired function must not resolve on the Spinach path.');
end

% Locate tracked sources independently of the working directory
this_file=[mfilename('fullpath') '.m'];
root=fileparts(fileparts(fileparts(this_file)));
[status,listing]=system(sprintf('git -C "%s" ls-files -- "*.m"',root));
assert(status==0,'The tracked-source inventory could not be read.');
files=strsplit(strtrim(listing),newline); offenders={};
names=['(?:' strjoin(retired,'|') ')'];
pattern=['(?<![\w.])' names '\s*\(|@\s*' names '(?!\w)'];
for n=1:numel(files)
    filename=fullfile(root,files{n});
    if strcmp(filename,this_file), continue; end
    source=strip_text(fileread(filename));
    source=regexprep(source,'\.\.\.[^\r\n]*','');
    if ~isempty(regexp(source,pattern,'once'))
        offenders{end+1}=files{n}; %#ok<AGROW>
    end
end
result=test_true(result,'No executable legacy references',isempty(offenders),...
                 ['Tracked sources containing retired calls: ' strjoin(offenders,', ')]);

% Check the native calculation and experiment entry points used by examples
native={'clexpn','spin_bath','context_task','orientation_average','powder',...
        'free_evolution','hahn_echo','cpmg_echo','cce_eseem',...
        'bath_pumped_hahn','finite_drive','noise_correlation','refocusing_train',...
        'cce_spectrum'};
for n=1:numel(native)
    result=test_true(result,['Native ' native{n}],...
                     exist(native{n},'file')==2&&startsWith(which(native{n}),[root filesep]),...
                     'The native entry point must resolve inside this Spinach checkout.');
end

% Exercise comments, escaped quotes, strings, and transpose-adjacent calls
quoted=sprintf('%% cce(a)\nx=''cce(a)''; y="cce_bath(a)"; z=''it''''s cce(a)'';');
result=test_true(result,'Quoted references are inert',...
                 isempty(regexp(strip_text(quoted),pattern,'once')),...
                 'Comment and string contents are not executable references.');
transposed='x=A''; y=cce(a);';
result=test_true(result,'Transpose does not hide a call',...
                 ~isempty(regexp(strip_text(transposed),pattern,'once')),...
                 'An executable call after a transpose must still be found.');

end

% Remove comments and quoted text while retaining transpose operators
function source=strip_text(source)

% Remove block comments before scanning line comments and quoted literals
source=regexprep(source,'(?m)^\s*%\{[^\r\n]*[\s\S]*?^\s*%\}[^\r\n]*','');
n=1;
while n<=numel(source)

    % A line comment ends at the next line boundary
    if source(n)=='%'
        last=find(source(n:end)==newline,1);
        if isempty(last), last=numel(source)-n+2; end
        source(n:n+last-2)=' '; n=n+last-1;
        continue
    end

    % A quote after an adjacent expression token denotes transposition
    quote=source(n);
    if quote==char(39)&&n>1&&...
       ~isempty(regexp(source(n-1),'[\w)\]}.''"]','once'))
        n=n+1; continue
    end
    if quote==char(39)||quote=='"'
        first=n; n=n+1;
        while n<=numel(source)
            if source(n)==quote
                if n<numel(source)&&source(n+1)==quote
                    n=n+2; continue
                end
                n=n+1; break
            end
            n=n+1;
        end
        source(first:n-1)=' '; source(n-1)='0';
    else
        n=n+1;
    end

end

end


