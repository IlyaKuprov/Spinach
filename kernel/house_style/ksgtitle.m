% House style settings for Matlab figures; a product of much
% experience with academic publication aesthetics. Syntax:
%
%                         ksgtitle(x)
%
% Parameters:
%
%    x - a character string
%
% Outputs:
% 
%    creates or updates the overall title
%    of the current grid of plots
%
% a.acharya@soton.ac.uk
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=ksgtitle.m>

function ksgtitle(x)

% Check consistency
grumble(x);

% Bold title rendered by LaTeX
sgtitle(['\textbf{' x '}'],'Interpreter','latex');

end

% Consistency enforcement
function grumble(x)
if ~ischar(x)
    error('x must be a character string');
end
end

% "There is considerable overlap between the intelligence 
%  of the smartest bears and the dumbest tourists."
%
% A forest ranger at the Yosemite National 
% Park on why it is hard to design the per-
% fect grabage bin to keep the bears from
% breaking into it.

