% Docs header here.

function [alpha,fx]=cubic_interp(End_A,End_B,alpha_A,alpha_B,f_A,dir_deriv_A,f_B,dir_deriv_B)

% Get the coefficients of the cubic polynomial
c1=-2*(f_B-f_A)+(dir_deriv_A+dir_deriv_B)*(alpha_B-alpha_A);
c2= 3*(f_B-f_A)-(2*dir_deriv_A+dir_deriv_B)*(alpha_B-alpha_A);
c3=(alpha_B-alpha_A)*dir_deriv_A; c4=f_A;

% Convert bounds to the z-space
bounds=([End_A End_B]-alpha_A)/(alpha_B-alpha_A);

% Find derivative roots
sPoints = roots([3*c1 2*c2 1*c3]);

% Remove complex roots
sPoints(imag(sPoints)~=0)=[];

% Remove points outside ranges 
sPoints(sPoints<min(bounds))=[];
sPoints(sPoints>max(bounds))=[];

% Make vector with solutions
sPoints=[min(bounds) sPoints(:)' max(bounds)];

% Select the global maximum point
[fx,k]=max(polyval([c1 c2 c3 c4],sPoints));

% Add the offset and scale back from [0..1] to the alpha domain
alpha = alpha_A + sPoints(k)*(alpha_B - alpha_A);

end

% It is baked beans that should be banned.
% In nauseating sauce encanned,
% All sickly-sweet and liquidous
% The bastards are ubiquitous.
%
% An English Breakfast should eschew
% This transatlantic parvenu –
% No beans! No beans! I cry. Too late.
% Plop goes the dollop on my plate.
%
% Not even a mixed grill avoids
% The devil's horrid haemorrhoids.
% I scrape them to the side, but still
% They manage to pollute my meal.
%
% They sidle up and touch my chips,
% Contriving thus to reach my lips
% So I again am forced to taste
% Their overtones of toxic waste.
%
% Ann Drysdale

