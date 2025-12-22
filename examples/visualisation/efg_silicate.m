% Example of electric field gradient tensor visualisation for
% an aluminosilicate solid. CASTEP log is parsed.
%
% ilya.kuprov@weizmann.ac.il
% e.dib@soton.ac.uk
% m.carravetta@soton.ac.uk

function efg_silicate()

% Import CASTEP data
props=c2spinach('alsilicate.magres');

% Do the visualisation
kfigure(); scale_figure([1.875 1.125]);
subplot(1,2,1); options.style='ellipsoids';
efg_display(props,{'Al'},100,[],options);
set(gca,'CameraPosition',[40 40 40]);
ktitle('ellipsoids');
subplot(1,2,2); options.style='harmonics';
efg_display(props,{'Al'},100,[],options);
set(gca,'CameraPosition',[40 40 40]);
ktitle('spherical harmonics');

end

