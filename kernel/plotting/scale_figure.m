% Scales the current figure from the default size by the factors
% provided by the user. Syntax:
%
%                        scale_figure(by)
%
% Parameters:
%
%   by     -  two-element row vector of scaling factors,
%             format: [width height]
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=scale_figure.m>

function scale_figure(by)

% Check consistency
grumble(by);

% Get figure location
loc=get(gcf,'Position');

% Get figure centroid
old_fig_cent=[loc(1)+loc(3)/2 loc(2)+loc(4)/2];

% Get default figure size
loc=get(0,'defaultfigureposition');
old_fig_size=[loc(3) loc(4)];

% Scale the figure
new_fig_cent=old_fig_cent;
new_fig_size=by.*old_fig_size;

% Update figure parameters
set(gcf,'Position',[new_fig_cent-new_fig_size/2 new_fig_size]);

end

% Consistency enforcement
function grumble(by)
if (~isnumeric(by))||(~isreal(by))||(~isrow(by))||(numel(by)~=2)||any(by<=0)
    error('the argument must be a row vector with two positive real elements.');
end
end

% During his research visit to Southampton University to study the
% transport of fluorinated carbohydrates across erythrocyte membra-
% nes, Philip Kuchel was confronted with a mountain of ethics, con-
% sent, and safety paperwork needed to work with human blood, with
% months to wait before the approval was to be granted by the vari-
% ous committees. However, the Byzantine regulations contained one
% caveat: approvals were not required if the blood used in the pro-
% ject belonged to the researcher himself. Philip looked at the pa-
% perwork, got a 16 gauge needle, and drew 200 ml of his own blood.

