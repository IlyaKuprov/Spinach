% Example of carbon hyperfine tensor visualisation for
% pyrene cation radical. Gaussian log is parsed.
%
% ilya.kuprov@weizmann.ac.il

function hfc_pyrene()
    
% Read the Gaussian log
props=gparse('pyrene_cation.log');

% Do the visualization
kfigure(); 
subplot(1,2,1);
options.style='ellipsoids';
hfc_display(props,{'C'},0.2,[],options);
set(gca,'CameraPosition',[40 40 40]);
ktitle('ellipsoids');
subplot(1,2,2);
options.style='harmonics';
hfc_display(props,{'C'},0.2,[],options);
set(gca,'CameraPosition',[40 40 40]);
ktitle('spherical harmonics');
scale_figure([1.875 1.125]);

end

