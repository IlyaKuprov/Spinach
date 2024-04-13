% Extracting the susceptibility tensor from DFT hyperfine tensors and
% experimental paramagnetic shifts. A combinatorial procedure is used 
% that cycles through ambiguous assignments.
%
% Calculation time: minutes.
%
% i.kuprov@soton.ac.uk
% s.g.worswick@soton.ac.uk

function combi_fit_2()

% Read DFT HFCs in Gauss
props=gparse('l2_parker_funk_yb.log');
parameters.hfcs=props.hfc.full.matrix;
parameters.nel=1;

% Isotope list
parameters.isotopes=cell(27,1);
for n=1:27
    parameters.isotopes{n}='1H';
end

% Spin groups with identical PCS
parameters.spin_groups={[28 22 30]; [25 33 81]; [80 32 24]; [27 82 21]; 
                        [36 58 85]; [35 57 84]; [ 2 88 39]; [ 4 90 41]; [61 92 43]};   

% Diamagnetic shifts
parameters.d_shifts=[3.62 2.65 2.65 2.86 4.95 4.10 7.40 8.00 7.80];

% Diamagnetic shift ambiguities
parameters.d_ambig={[1 2 3 4],[5 6],[7 8 9]};

% Paramagnetic shifts
parameters.p_shifts=[3.8  -13.7   -0.6   -3.4    5.2   20.7   10.7   10.3   10.7];

% Paramagnetic shift ambiguities
parameters.p_ambig={[3 4]};

% Run the combinatorial fitting
[~,~,pcs_theo,pcs_expt,~,~]=pcs_combi_fit(parameters);

% Plot theory against experiment
figure();
plot(pcs_theo,pcs_expt,'ro'); hold on;
plot([-20 20],[-20 20],'b-'); hold off;
xlabel('Theoretical PCS / ppm');
ylabel('Experimental PCS / ppm');
axis tight; kgrid;

end

