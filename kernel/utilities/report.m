% Writes a log message to the console or an ACSII file. The message
% includes the call stack of the function that produced it. Syntax:
%
%                report(spin_system,report_string)
%
% Parameters:
%
%    report_string  - a character string 
%
% Outputs:
%
%    this function prints the message to the console or to the
%    destination specified in spin_system.sys.output
%
% Note: a newline symbol at the end of the string is not neces-
%       sary - it is added by the function.
%
% Note: all output produced by this function may be silenced 
%       by setting sys.output='hush' in the Spinach input 
%       stream or by setting spin_system.sys.output='hush'
%       at any point during the calculation.
%
% i.kuprov@soton.ac.uk
% ledwards@cbs.mpg.de
%
% <https://spindynamics.org/wiki/index.php?title=report.m>

function report(spin_system,report_string)

% Catch single-argument calls
if nargin==1
    error('console reporting function requires two arguments.');
end

% Ignore the call if the system is hushed
if ~strcmp(spin_system.sys.output,'hush')
    
    % Validate the input
    grumble(spin_system,report_string);
    
    % List call stack exceptions
    not_useful={'make_general_channel/channel_general','parallel_function'};
    mdcs={'remoteParallelFunction','remoteBlockExecution','cppRemoteParallelFunction'};
    
    % Get the call stack
    call_stack=dbstack;
    
    % Compose the prefix
    for n=1:numel(call_stack)

        % Uninformative entries
        if ismember(call_stack(n).name,not_useful)
            
            % Simply delete
            call_stack(n).name='';
        
        % Parallelisation rigging entries
        elseif ismember(call_stack(n).name,mdcs)

            % Tell the user it's a parallel call
            call_stack(n).name=' (parfor/spmd) > ';

        else

            % Just add the entry
            call_stack(n).name=[call_stack(n).name ' > '];

        end

    end

    % Concatenate the names
    prefix_string=[call_stack(end:-1:2).name];

    % Remove m-file extensions
    prefix_string=prefix_string(1:(end-3));

    % Fix occasional empty prefixes
    if isempty(prefix_string), prefix_string=' '; end
    
    % Roll the prefix
    if numel(prefix_string)<50
        prefix_string=pad(prefix_string,50);
    else
        prefix_string=['...' prefix_string((end-46):end)];
    end
    
    % Add prefix to the report string
    report_string=['[' prefix_string ' ]  ' report_string];
    
    % Send the report string to the output, ignoring impossible writes
    try fprintf(spin_system.sys.output,'%s\n',report_string); end %#ok<TRYNC>
    
end

end

% Consistency enforcement
function grumble(spin_system,report_string)
if (~isfield(spin_system,'sys'))||~isfield(spin_system.sys,'output')
    error('spin_system.sys.output field must exist.');
end
if ((~isa(spin_system.sys.output,'double'))&&(~isa(spin_system.sys.output,'char')))||...
   (isa(spin_system.sys.output,'char')&&(~strcmp(spin_system.sys.output,'hush')))
    error('spin_system.sys.output must be either ''hush'' or a file ID.');
end
if ~ischar(report_string)
    error('report_string must be a string.');
end
end

% All parts should go together without forcing. You must remember that the
% parts you are reassembling were disassembled by you. Therefore, if you
% can't get them together again, there must be a reason. By all means, do
% not use a hammer.
%
% IBM Manual, 1925 

