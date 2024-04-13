% Interface to the Bruker SpinJet AWG, calling a library of Python
% scripts that inteface with the Bruker Xepr API. Syntax:
%
%           data=awg_interface(spin_system,awg_cmd,cmd_input)
%
% Parameters:
%
%       spin_system     - should contain spinach root directory
%                         in spin_system.sys.root_dir
%
%       awg_cmd         - switch between available commands:
%                               'reset_pspel'
%                               'compile_pspel_shp'
%                               'compile_pspel_def'
%                               'compile_pspel_exp'
%                               'modify_pspel_defs'
%                               'acquire_data'
%
%       cmd_input       - cell array of inputs to be passed to 
%                         python scripts
%
% Outputs:
%
%       data.X          - abscissa values from the acquired 
%                         signal
%
%       data.rY         - real ordinate values from the ac-
%                         quired signal
%
%       data.iY         - imaginary ordinate values from 
%                         the acquired signal.
%
% david.goodwin@inano.au.dk
%
% <https://spindynamics.org/wiki/index.php?title=awg_interface.m>

function data=awg_interface(spin_system,awg_cmd,cmd_input)

% Bootstrap cmd_input if not supplied
if ~exist('cmd_input','var'), cmd_input={}; end

% Check consistency
grumble(spin_system,awg_cmd,cmd_input);

% Initialise data structure
data=[];

% Switch between defined routines
switch awg_cmd
    
    case {'reset_pspel'}
        
        % Reset pulse-spel by connecting with Xepr
        py_run(spin_system,'Xepr_resetexpt');
        
        % Instruct user to hide then show the PulseSPEL window (Bruker bug)
        report(spin_system,'<<user: 1) button press hide PulseSPEL window>>');
        report(spin_system,'<<user: 2) button press show PulseSPEL window>>');
        report(spin_system,'...pausing... keypress to continue (or "ctrl+c" to end simulation)...');
        pause
        
    case {'compile_pspel_shp'}
        
        % Compile the shape file defined in cmd_input
        py_run(spin_system,'Xepr_plsspel_shpfile',cmd_input);
        
    case {'compile_pspel_def'}
        
        % Compile the definitions file defined in cmd_input
        py_run(spin_system,'Xepr_plsspel_deffile',cmd_input);
        
    case {'ompile_pspel_exp'}
        
        % Compile the experiment file defined in cmd_input
        py_run(spin_system,'Xepr_plsspel_expfile',cmd_input);
        
    case {'modify_pspel_defs'}
        
        % Run command to modify definitions in definitions file
        py_run(spin_system,'Xepr_plsspel_moddefs',cmd_input);
        
    case {'acquire_data'}
        
        % Acquire data from spectometer
        py_run(spin_system,'Xepr_getdata',cmd_input);
        
        % Save x axis, imaginary, and real parts of the signal
        data.X  = dlmread([cmd_input{3} '_X.txt']);  %#ok<DLMRD> 
        data.rY = dlmread([cmd_input{3} '_rY.txt']); %#ok<DLMRD> 
        data.iY = dlmread([cmd_input{3} '_iY.txt']); %#ok<DLMRD> 
        
        % Delete temporary file of results
        delete([cmd_input{3} '_X.txt'],[cmd_input{3} '_rY.txt'],[cmd_input{3} '_iY.txt'])
        
end

end

% Consistency enforcement
function grumble(spin_system,awg_cmd,cmd_input)
if ~isfield(spin_system,'sys') || ~isfield(spin_system.sys,'root_dir')
    error('spin_system.sys.root_dir must be supplied')
elseif ~isfolder(spin_system.sys.root_dir)
    error('spin_system.sys.root_dir should be a path to a directory')
end
if ~ischar(awg_cmd) || ~ismember(awg_cmd,{'reset_pspel','compile_pspel_shp',...
        'compile_pspel_def','compile_pspel_exp','modify_pspel_defs','acquire_data'})
    error('unknown awg_cmd command')
end
if ~iscell(cmd_input)
    error('cmd_input should be supplied as a cell array of python script inputs')
end
if ismember(awg_cmd,{'acquire_data'}) && numel(cmd_input)~=3
    error('cmd_input should contain 3 strings when acquire_data is used')
end
if ismember(awg_cmd,{'modify_pspel_defs'}) && (numel(cmd_input)<2 || mod(numel(cmd_inputs),2)~=0)
    error('must provide input cell array of even number of elements to modify definitions')
end
if ismember(awg_cmd,{'compile_pspel_exp'}) && numel(cmd_input)~=1
    error('exp_file should be suppied in cmd_input.')
end
if ismember(awg_cmd,{'compile_pspel_def'}) && numel(cmd_input)~=1
    error('def_file should be suppied in cmd_input.')
end
if ismember(awg_cmd,{'compile_pspel_shp'}) && numel(cmd_input)~=1
    error('shp_file should be suppied in cmd_input.')
elseif ismember(awg_cmd,{'compile_pspel_shp'})
    shp_file_info=dir(cmd_input{1});
    if shp_file_info.bytes>=262144
        error(['shape file too large on disk (' int2str(shp_file_info.bytes) ' bytes) - maximum 262144 bytes'])
    end
end
end

% If you want to help a novice - work with him.
% If you want to help an old man - work for him.
% If you want to help a master - don't get in the way.
% If you want to help a fool - you are a fool too.
%
% A Russian saying

