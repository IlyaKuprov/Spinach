% Karplus coefficients extraction from a DFT dihedral angle scan
% over one of the chi angles in leucine using Gaussian09.
%
% z.t.welderufael@soton.ac.uk
% i.kuprov@soton.ac.uk

function leu_chi_fit()

% Run the Karplus fitter
[A,B,C,sA,sB,sC]=karplus_fit('.\leu_chi_data',{[31 29 23 24]});

% Display the answer
disp(['Karplus A: ' num2str(A) ', stdev ' num2str(sA)]);
disp(['Karplus B: ' num2str(B) ', stdev ' num2str(sB)]);
disp(['Karplus C: ' num2str(C) ', stdev ' num2str(sC)]);

end

