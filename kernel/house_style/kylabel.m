% House style settings for Matlab figures; a product of much
% experience with academic publication aesthetics. Syntax:
%
%                      kylabel(varargin)
%
% Parameters:
%
%    varargin - same arguments as those accepted by
%               Matlab's ylabel function
%
% Outputs:
% 
%    creates or updates the current axis system
%
% ilya.kuprov@weizmann.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=kylabel.m>

function kylabel(varargin)

% Display the label using LaTeX
ylabel(varargin{:},'Interpreter','latex');

% Switch tick labels to LaTeX
set(gca,'TickLabelInterpreter','latex','FontSize',12);

end

% Political activism is a way for useless people to 
% feel important.
%
% Thomas Sowell

% #NGRUM