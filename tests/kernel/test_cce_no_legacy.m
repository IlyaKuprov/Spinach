% Absence of the retired CCE interface and its executable callers. Syntax:
%
%                    result=test_cce_no_legacy()
%
% Outputs:
%
%     result - scalar Spinach regression-test result structure
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
         'cce_observe','cce_prepare','cce_sequence','cce_validate','cce_controls',...
         'context_task','orientation_average'};
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
files=strsplit(strtrim(listing),newline);
offenders={};
names=['(?:' strjoin(retired,'|') ')'];
pattern=['(?<![\w.])' names '\s*\(|@\s*' names '(?!\w)'];
for n=1:numel(files)
    filename=fullfile(root,files{n});
    if strcmp(filename,this_file)
        continue;
    end
    source=strip_text(fileread(filename));
    source=regexprep(source,'\.\.\.[^\r\n]*','');
    if ~isempty(regexp(source,pattern,'once'))
        offenders{end+1}=files{n}; %#ok<AGROW>
    end
end
result=test_true(result,'No executable legacy references',isempty(offenders),...
                 ['Tracked sources containing retired calls: ' strjoin(offenders,', ')]);

% Check the native calculation and experiment entry points used by examples
native={'clexpn','spin_bath','powder',...
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

% Every retired parameter path is rejected by the native input boundary
sys.isotopes={'1H','13C'};
sys.magnet=0;
sys.output='hush';
sys.parallel={'processes',1};
sys.disable={'hygiene'};
inter.zeeman.matrix={zeros(3),zeros(3)};
spin_system=create(sys,inter);
p.observed=1;
p.times=[0 0.1];
p.transition=struct('spins',1,'ms',[0.5 -0.5]);
p.cce=struct('approximation','generalized','order',1);
p.cce.numerics=struct('atol',1e-12,'rtol',1e-10,'max_memory_bytes',2^28);
paths={'bath_state.temperature',...
       'cluster_cap',...
       'clusters',...
       'controls',...
       'deer.enabled',...
       'diagnostics',...
       'initial_rho',...
       'method',...
       'numerics.full_bath',...
       'quantity',...
       'run',...
       'sampling',...
       'sequence',...
       'sequence.amplitude',...
       'sequence.angle',...
       'sequence.at',...
       'sequence.axis',...
       'sequence.count',...
       'sequence.detuning',...
       'sequence.duration',...
       'sequence.name',...
       'sequence.operator',...
       'sequence.placement',...
       'sequence.segments',...
       'sequence.target',...
       'sequence.type',...
       'cce.bath_state.temperature',...
       'cce.numerics.full_bath',...
       'cce.cluster_cap',...
       'cce.clusters',...
       'pump.enabled'};
for n=1:numel(paths)
    names=strsplit(paths{n},'.');
    request=p;
    request=setfield(request,names{:},1);
    rejected=false;
    try
        clexpn(spin_system,@free_evolution,request,'labframe');
    catch exception
        rejected=startsWith(exception.identifier,'cce_options:');
    end
    result=test_true(result,['Removed input ' paths{n}],rejected,...
                     'Removed fields must fail before constructing a retained generator.');
end

% Direct native experiments also reject raw controls, including empty arrays
[p,~]=cce_options(spin_system,p,@free_evolution,'labframe');
model=cce_compile(spin_system,p,'labframe');
[cluster,~]=cce_subsystem(model,2,struct(),[]);
[H,R,K,p]=cce_generators(cluster,p,[]);
p.nloops=1;
p.drive=struct('spins',1,'axis','x','amplitude',1,'detuning',0,'duration',0.1);
p.refocusing={struct('spins',1,'axis','x','angle',pi,'fraction',0.5)};
experiments={@free_evolution,@finite_drive,@refocusing_train,@hahn_echo,...
             @cpmg_echo,@cce_eseem,@bath_pumped_hahn,@noise_correlation};
for n=1:numel(experiments)
    for value={[],struct([]),struct('type','finite'),struct('type','ideal')}
        request=p;
        request.sequence=value{1};
        rejected=false;
        try
            experiments{n}(cluster,request,H,R,K);
        catch exception
            rejected=contains(exception.message,'finite_drive')||...
                     contains(exception.message,'refocusing_train');
        end
        result=test_true(result,['Raw events rejected by ' func2str(experiments{n})],...
                         rejected,'The native experiment must name the supported control interface.');
    end
end

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
        if isempty(last)
            last=numel(source)-n+2;
        end
        source(n:n+last-2)=' ';
        n=n+last-1;
        continue;
    end

    % A quote after an adjacent expression token denotes transposition
    quote=source(n);
    if quote==char(39)&&n>1&&...
       ~isempty(regexp(source(n-1),'[\w)\]}.''"]','once'))
        n=n+1;
        continue;
    end
    if quote==char(39)||quote=='"'
        first=n;
        n=n+1;
        while n<=numel(source)
            if source(n)==quote
                if n<numel(source)&&source(n+1)==quote
                    n=n+2;
                    continue;
                end
                n=n+1;
                break;
            end
            n=n+1;
        end
        source(first:n-1)=' ';
        source(n-1)='0';
    else
        n=n+1;
    end

end

end


