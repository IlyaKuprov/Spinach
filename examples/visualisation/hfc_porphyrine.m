% Example of proton hyperfine tensor visualisation for
% copper porphyrine. ORCA log is parsed.
%
% ilya.kuprov@weizmann.ac.il
% e.suturina@soton.ac.uk

function hfc_porphyrine()
    
% Read the Gaussian log
props=oparse('porphyrine.out');

% Do the visualisation
figure(); subplot(1,2,1);
options.style='ellipsoids';
hfc_display(props,{'H'},2.0,[],options);
set(gca,'CameraPosition',[40 40 40]);
ktitle('ellipsoids');
subplot(1,2,2);
options.style='harmonics';
hfc_display(props,{'H'},2.0,[],options');
set(gca,'CameraPosition',[40 40 40]);
ktitle('spherical harmonics');
scale_figure([1.875 1.125]);

end

