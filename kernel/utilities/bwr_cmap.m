% Blue -> White -> Red colour map with 255 points and 
% white colour corresponding to zero. Syntax:
%
%                      cmap=bwr_cmap()
%
% The output is 255x3 RGB column that starts at blue,
% goes into white and then into red in with a quadra-
% tic bend.
%
% Outputs:
%
%    cmap - colour map in Matlab format
%   
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=bwr_cmap.m>

function cmap=bwr_cmap()

% Preallocate the map
cmap=zeros(255,3);

% Rise from blue to white
cmap(1:128,1)=linspace(0,1,128); % reds
cmap(1:128,2)=linspace(0,1,128); % greens
cmap(1:128,3)=1;                 % blues

% Rise from white to red
cmap(128:255,1)=1;                         % reds
cmap(128:255,2)=fliplr(linspace(0,1,128)); % greens
cmap(128:255,3)=fliplr(linspace(0,1,128)); % blues

% Improve contrast
cmap=cmap.^2;
    
end

% The worst thing I can be is the same as everybody 
% else. I hate that. 
%
% Arnold Schwarzenegger

