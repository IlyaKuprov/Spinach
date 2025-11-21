% Example of carbon shielding tensor visualisation for
% strychnine molecule. Gaussian log is parsed.
%
% Note: antisymmetric components of the shielding ten-
%       sors are ignored.
%
% ilya.kuprov@weizmann.ac.il

function cst_strychnine()
    
% Read the Gaussian log
props=gparse('strychnine.log');

% Do the visualisation
kfigure(); subplot(1,2,1);
options.style='ellipsoids';
cst_display(props,{'C'},0.005,[],options);
set(gca,'CameraPosition',[40 40 40]);
ktitle('ellipsoids');
subplot(1,2,2);
options.style='harmonics';
cst_display(props,{'C'},0.01,[],options);
set(gca,'CameraPosition',[40 40 40]);
ktitle('spherical harmonics');
scale_figure([1.875 1.125]);

end

