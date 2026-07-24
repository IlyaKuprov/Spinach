% Field dependence of proton relaxation near the dysprosium centre
% in the DyL1 complex. Molecular and relaxation parameters correspond
% to a room-temperature proton NMR experiment.
%
% Calculation time: seconds
%
% ilya.kuprov@weizmann.ac.il

function sbm_dyl1_h1()

% Experiment parameters
fields=1.0:0.25:20.0;
nucleus='1H';

% Paramagnetic centre parameters
e_spin=15/2;
g_eff=4/3;
tau_r=140e-12;
tau_e=[0.2e-12 0.2e-12];

% H1 position relative to the metal centre
dist=4.3702201672;
a_iso=0;

% Combined correlation times
tau_c=1./(1/tau_r+1./tau_e);

% Preallocate relaxation rates
r1=zeros(size(fields));
r2=zeros(size(fields));

% Compute field dependence
for n=1:numel(fields)
    [r1_parts,r2_parts]=rlx_sbm(fields(n),nucleus,dist,a_iso,...
                                e_spin,g_eff,tau_c,tau_e);
    r1(n)=sum(r1_parts);
    r2(n)=sum(r2_parts);
end

% Plot relaxation rates
kfigure();
semilogy(fields,r1,'b-',fields,r2,'r-');
kxlabel('Magnetic field, Tesla');
kylabel('Relaxation rate, Hz');
klegend({'R_1','R_2'},'Location','northwest');
kgrid();

end


