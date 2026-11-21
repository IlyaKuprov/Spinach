% Two-frequency MAXCRAZED benchmark for water (W) and acetone (A).
% Chen, Branca and Warren, J. Chem. Phys. 136, 204509 (2012),
% doi:10.1063/1.4721637, Fig. 1(c) and Table I, d2=60 ms.
%
% The source does not give an absolute gradient calibration, mixture
% diffusivities, or spectral processing for Table I. This unfitted
% prediction uses a periodic bulk field, a nominal 1 mm X grating,
% reflecting diffusion on the simulation box, and no relaxation or
% radiation damping. See the accompanying knowledge entry for details.
%
% The result structure contains the complex FIDs and magnitude spectra
% (indirect frequency x direct frequency x X/XY/XZ), frequency axes (Hz),
% peak comparison table, and measured and predicted hetero/homo ratios.

function result=ddf_maxcrazed()

% Compare the nominal 60 ms inter-lobe gap with the Table I experiment
result=ddf_maxcrazed_scan([12 16 16],0.5e-3,60e-3,0);
f1_peaks=result.f1_peaks;
f2_peaks=result.f2_peaks;
predicted=result.predicted;
measured=[1.00 0.23 0.07;...
          0.44 0.11 0.12;...
          0.51 0.15 0.18;...
          0.48 0.14 0.04];

% Report residuals before ratios so the absolute suppression remains visible
result.peaks=table(f1_peaks,f2_peaks,measured,predicted,predicted-measured,...
                   'VariableNames',{'f1_hz','f2_hz','measured','predicted','residual'},...
                   'RowNames',{'2W,W','A+W,W','A+W,A','2A,A'});
result.measured_ratio=[0.64 0.68 2.65];
result.predicted_ratio=sum(predicted(2:3,:),1)./sum(predicted([1 4],:),1);
disp('Table I comparison: matrix columns are X, XY, XZ');
disp(result.peaks);
disp('Hetero/homo ratios: measurement, prediction, and residual');
disp([result.measured_ratio;result.predicted_ratio;...
      result.predicted_ratio-result.measured_ratio]);

end


