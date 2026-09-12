% Draws an academic journal style letter label
% in the top left corner of the current axis
% set. Syntax:
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
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=pauli.m>

function kletter(letter_label)

% Just place the label
text(0.025,0.875,letter_label,'Units','normalized',...
     'FontWeight','bold','FontSize',16);

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

