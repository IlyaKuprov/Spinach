% A wrapper intended to trick SPMD blocks into saving data. Can
% only save one variable at a time, its name in the mat file is
% "variable". Syntax:
%
%                 save_anyway(file_name,variable)
%
% Parameters:
%
%    file_name  - a character string specifying the 
%                 file name
%
%    variable   - the variable to be saved
%
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=save_anyway.m>

function save_anyway(file_name,variable)

% Check consistency
grumble(file_name);

% Just call save
save(file_name,'variable','-v7.3');

end

% Consistncy enforcement
function grumble(file_name)
if ~ischar(file_name)
    error('file_name must be a character string.');
end
end

% Life struggles to survive here, and while some clings 
% to a tenacious existence, it is anemic and sickly.
%
% Stellaris - a computer game

