% House style settings for Matlab figures; a product of much
% experience with academic publication aesthetics. Syntax:
%
%                      kzlabel(varargin)
%
% Parameters:
%
%    varargin - same arguments as those accepted by
%               Matlab's zlabel function
%
% Outputs:
% 
%    creates or updates the current axis system
%
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=kzlabel.m>

function kzlabel(varargin)

% Display the label using LaTeX
zlabel(varargin{:},'Interpreter','latex');

% Switch tick labels to LaTeX
set(gca,'TickLabelInterpreter','latex','FontSize',12);

end

% In an age of crybabies and professional victims,
% Rupert stood out like a saint in hell.
%
% Taki Theodoracopulos, 
% about Rupert Hambro

% #NGRUM