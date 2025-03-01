% House style settings for Matlab figures; a product of much
% experience with academic publication aesthetics. Syntax:
%
%                      klegend(varargin)
%
% Parameters:
%
%    varargin - same arguments as those accepted by
%               Matlab's legend function
%
% Outputs:
% 
%    creates or updates the current figure legend
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=klegend.m>

function klegend(varargin)

% Display the legend using LaTeX
leg_obj=legend(varargin{:},'Interpreter','latex');

% Make legend box translucent
set(leg_obj.BoxFace,'ColorType','truecoloralpha',...
                    'ColorData',uint8([200 200 200 64]'));

end

% It was awesome - my first tabloid story. If you're going to
% have a tabloid story written about you, it might as well be
% with Johnny Depp.
%
% Christina Ricci, about newspapers 
% reporting her romance with Johnny 
% Depp

% #NGRUM