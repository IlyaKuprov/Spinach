% Draws an academic journal style letter label in the
% top left corner of the current axis set. The label is
% placed inside the outer box of the axes (the whole fi-
% gure when there is one axis set), with its left edge
% and its cap line 10 points away from the left and the
% top edge of that box, whatever the proportions of the
% figure. Syntax:
%
%             kletter(letter_label)
%
% Parameters:
%
%   letter_label - one-element character
%                  string with the label
%
% Outputs:
%
%   updates the current axis system
%
% Note: the offsets are measured when the function is
%       called; resize the figure before labelling it.
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=kletter.m>

function kletter(letter_label)

% Check consistency
grumble(letter_label);

% Offset from the box edges, points
edge_offset=10;

% Plot box and outer box of the current axes in points
ax_obj=gca; ax_units=ax_obj.Units; ax_obj.Units='points';
plot_box=ax_obj.Position; outer_box=ax_obj.OuterPosition;
ax_obj.Units=ax_units;

% Label position relative to the plot box corner
label_x=outer_box(1)-plot_box(1)+edge_offset;
label_y=outer_box(2)+outer_box(4)-plot_box(2)-edge_offset;

% Place the label with its cap line at the top offset
text(label_x,label_y,letter_label,'Units','points',...
     'HorizontalAlignment','left','VerticalAlignment','cap',...
     'FontWeight','bold','FontSize',16);

end

% Consistency enforcement
function grumble(letter_label)
if (~ischar(letter_label))||(~isscalar(letter_label))
    error('letter_label must be a one-element character string.');
end
end

%                                    .-"""-.
%                                   / .--.  \
%                                  | |    \_/
%                                   \ \
%                                     \ \
%      /\           /\                 | |
%     /  \_________/  \    _____       / /
%    /   ,         ,   \.-'     `-.___/ /
%   /    o         o    \              /
%   |         Y         |              \
%    \      \___/      /                |
%     `-.__       __.-'                 |
%          /`-----'          ___       /
%         /      |         .'   `.    /
%         |  |   |________/       \  /
%         |  |   |         |      / /
%         |  |   |         |     / /
%         |  |   |         |    / /
%        (___|___)        (____|_)
%    +---------------------------------------+
%    |                                       |
%    |        По техническим причинам        |
%    |       кот нассал вам в капучино       |
%    |                                       |
%    +---------------------------------------+

