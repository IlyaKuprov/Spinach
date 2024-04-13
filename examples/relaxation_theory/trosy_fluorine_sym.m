% Transverse relaxation rate as a function of the applied magnetic 
% field in a 3-fluorotyrosine labelled protein. The fluorine atom
% and its directly bonded carbon are included. Analytical calcula-
% tions broken down by mechanism.
%
% Calculation time: seconds.
% 
% i.kuprov@soton.ac.uk

function trosy_fluorine_sym()

% Read 3-fluorotyrosine DFT calculation
[~,inter_dft]=g2spinach(gparse('../standard_systems/3_fluoro_tyr.log'),...
                            {{'C','13C'},{'F','19F'}},[186.38 192.97]);
                                 
% Extract coordinates and CSAs
ZF=inter_dft.zeeman.matrix{8};
ZC=inter_dft.zeeman.matrix{7};
RF=inter_dft.coordinates{8};
RC=inter_dft.coordinates{7};

% Magnetic field grid
lin_freq=linspace(200,800,20);
B0=2*pi*lin_freq*1e6/spin('1H');

% Loop over magnetic fields
for n=1:numel(B0) %#ok<*AGROW>
    
    % Call the analytical function
    [F,C]=rlx_dd_csa(B0(n),25e-9,{'19F','13C'},{ZF,ZC},{RF,RC});
      
    % Relaxation rates
    r2c(n)=C.r2.total; 
    r2f(n)=F.r2.total;
    f_bro(n)=F.trosy.total_bro;
    f_nar(n)=F.trosy.total_nar;
    c_bro(n)=C.trosy.total_bro;
    c_nar(n)=C.trosy.total_nar;
    
    % Mechanisms for 13C
    c_tro_dd(n)=C.trosy.dd;
    c_tro_csa(n)=C.trosy.csa;
    c_tro_xc(n)=C.trosy.xc;
    
end
                         
% Plotting
figure();
plot(lin_freq',[f_bro' r2f' f_nar']); ylim([0 1200]);
kxlabel('Proton Larmor frequency, MHz'); kgrid;
kylabel('Relaxation matrix element, Hz');
klegend({'TROSY, broad','$R_2$ relax. rate','TROSY, narrow'},...
         'Location','northwest','FontSize',12);
ktitle('3-fluoro-Tyr, $^{19}$F TROSY line relaxation rates');
    
figure();
plot(lin_freq',[c_bro' r2c' c_nar']); ylim([0 140]);
kxlabel('Proton Larmor frequency, MHz'); kgrid;
kylabel('Relaxation matrix element, Hz');
klegend({'TROSY, broad',' $R_2$ relax. rate',' TROSY, narrow'},...
         'Location','northwest','FontSize',12);
ktitle('3-fluoro-Tyr, $^{13}$C TROSY line relaxation rates');
    
% TROSY rate by mechanism
figure(); 
bar(lin_freq,[c_tro_dd' c_tro_csa' -abs(c_tro_xc')],'stacked'); 
kxlabel('Proton Larmor frequency, MHz'); kgrid;
kylabel('Relaxation matrix element, Hz');
klegend({'DD','CSA','DD-CSA'},'Location','northwest','FontSize',12);
ktitle('3-fluoro-Tyr, $^{13}$C TROSY rate by mechanism');

end

                         