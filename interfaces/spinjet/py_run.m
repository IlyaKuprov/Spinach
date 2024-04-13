% Runs a python script from /interfaces/spinjet/Xepr_python/ folder of
% Bruker Xepr installation. Variable inputs and output can be passed to
% and from the script. Syntax:
%
%              arg_out=py_run(spin_system,pyscript,arg_in)
%
% Parameters:
%
%       spin_system   - spin_system created by spinach, with
%                       spin_system.sys.root_dir defined
%
%       pyscript      - name of the python script to run, in 
%                       the form of a string without the .py
%                       extension
%                       
%       arg_in        - cell array containing the inputs to 
%                       the script
%
% Outputs:
%
%       arg_out       - cell array of outputs from the script,
%                       split at every space character (only
%                       if no error exception has occourred -
%                       in that case the error message would
%                       be returned)
%
% david.goodwin@inano.au.dk
%
% <https://spindynamics.org/wiki/index.php?title=py_run.m>

function arg_out=py_run(spin_system,pyscript,arg_in)

% Initialise arg_in if not an input
if ~exist('arg_in','var'), arg_in={}; end

% Check consistency
grumble(spin_system,pyscript,arg_in);

% Initialise output as an empty cell
arg_out={};

% Initialise command string to run the selected python script file
command_string=['python ' spin_system.sys.root_dir ...
                '/interfaces/spinjet/Xepr_python/' pyscript '.py'];

% Process extra argument cell, if needed
if ~isempty(arg_in)
    
    % Append command string with space
    command_string=[command_string, ' '];
    
    % Process each string in the cell of inputs
    for inputs=1:numel(arg_in)
        
        % Ensure input as characters
        if ~ischar(arg_in{inputs})
            arg_in{inputs}=num2str(arg_in{inputs}); 
        end
        
        % Append the command string
        command_string=[command_string arg_in{inputs} ' ']; %#ok<AGROW>
        
    end
    
end

% Run Xepr python script with system command-line
[status, output_string] = system(command_string);

% Detect error exception within python script - if no error, then compile
% outputs from python script, if any were assigned during python print.
if ~(status==0) 
    error(['Exception with error code (' int2str(status) ') in ' pyscript ...
           '.py:\npython script returned with error > ' output_string],class(status))
elseif ~isempty(output_string)
    arg_out=strsplit(strtrim(output_string))';
end

end

% Consistency enforcement
function grumble(spin_system,py_script_str,python_inputs)
if isempty(spin_system)||(~isfield(spin_system,'sys'))||...
   (~isfield(spin_system.sys,'root_dir'))||(~ischar(spin_system.sys.root_dir))
    error('spinach root directory must be supplied within spin_system.sys.root')
end
if isempty(py_script_str)||(~ischar(py_script_str))
   if ~exist([spin_system.sys.root_dir '/interfaces/spinjet/Xepr_python/' py_script_str '.py'],'file') 
       error([spin_system.sys.root_dir '/interfaces/spinjet/Xepr_python/' py_script_str '.py does not exist'])
   end
end
if ~iscell(python_inputs)
    error(['python_inputs should be a cell array of all inputs to python script ' py_script_str '.py'])
end
end

% In 1847, Robert Liston performed an amputation in 150 seconds,
% operating so quickly that he accidentally amputated his assis-
% tant's fingers as well. Both the patient and the assistant la-
% ter died of sepsis. He also slashed through the coat tails of
% a distinguished surgical spectator, who was so terrified that
% the knife had pierced his vitals he dropped dead from fright.
% That was the only operation in history with a 300 percent mor-
% tality rate.
%
% On a different day, Liston had an argument with his house-sur-
% geon. Was the red, pulsating tumour in a small boy's neck a 
% straightforward abscess of the skin, or a dangerous aneurism
% of the carotid artery? 'Pooh!' Liston exclaimed impatiently.
% 'Whoever heard of an aneurism in one so young?' Flashing a 
% knife from his waistcoat pocket, he lanced it. Out leaped ar-
% terial blood, and the boy fell. The patient died but the ar-
% tery lives, in University College Hospital pathology museum,
% specimen No. 1256.
%
% Richard Gordon

