% Example of shielding tensor visualisation for a peptide
% bond. Gaussian log is parsed.
%
% Note: antisymmetric components of the shielding ten-
%       sors are ignored.
%
% i.kuprov@soton.ac.uk

function cst_peptide_bond()
    
% Read the Gaussian log
props=gparse('..\standard_systems\amino_acids\ala.log');

% Do the visualisation
figure(); scale_figure([2.0 1.0]);

subplot(1,3,1); options.style='harmonics';
cst_display(props,{'C'},0.01,[],options);
set(gca,'CameraPosition',[40 40 40]);
ktitle('$^{13}$C CST');

subplot(1,3,2); options.style='harmonics';
cst_display(props,{'H'},0.05,[],options);
set(gca,'CameraPosition',[40 40 40]);
ktitle('$^{1}$H CST');

subplot(1,3,3); options.style='harmonics';
cst_display(props,{'N'},0.01,[],options);
set(gca,'CameraPosition',[40 40 40]);
ktitle('$^{15}$N CST');

end

