% House style settings for Matlab figures; a product of much
% experience with academic publication aesthetics. Syntax:
%
%                      kxlabel(varargin)
%
% Parameters:
%
%    varargin - same arguments as those accepted by
%               Matlab's xlabel function
%
% Outputs:
% 
%    creates or updates the current axis system
%
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=kxlabel.m>

function kxlabel(varargin)

% Display the label using LaTeX
xlabel(varargin{:},'Interpreter','latex');

% Switch tick labels to LaTeX
set(gca,'TickLabelInterpreter','latex','FontSize',12);

end

% МОСКВА, 14 мая 2021 - РИА Новости: Полиция задержала жителя
% подмосковных Химок, бросившего телевизор в Вечный огонь, он
% был в состоянии наркотического опьянения.

% #NGRUM