% Example Diels-Alder cycloaddition reaction settings: pentadiene 
% (reactant), acrylonitrile (reactant), exo-norbornene (product),
% endo-norbornene (product), and cyanomethane (solvent). Atom co-
% ordinates and chemical shift anisotropies are pulled from a DFT
% calculation; isotropic chemical shifts and J-couplings are ex-
% perimental (some J-coupling signs are missing).
%
% Parameters:
%
%     none
%
% Outputs:
%
%     sys, inter, bas - Spinach input data structures, remember
%                       to specify the field in sys.magnet
%
% a.acharya@soton.ac.uk
% bruno.linclau@ugent.be
% madhukar.said@ugent.be
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=dac_reaction.m>

function [sys,inter,bas]=dac_reaction()

% Find own location
own_folder=mfilename('fullpath');
own_folder=own_folder(1:(end-12));

% Cyclopentadiene (substance A)
props_a=gparse([own_folder 'cyclopentadiene.log']);
[sys_a,inter_a]=g2spinach(props_a,{{'H','1H'}},31.8);

% Replace isotropic shifts 
inter_a.zeeman.matrix=shift_iso(inter_a.zeeman.matrix,[1 2 3 4 5 6],...
                               [6.628 6.427 6.427 6.628 2.797 2.797]);

% Create labels 
sys_a.labels={'HB','HA','HAp','HBp','HCp','HC'};

% Replace J-couplings 
inter_a.coupling.scalar=zeros(6,6);
inter_a.coupling.scalar(idxof(sys_a,'HA'), idxof(sys_a,'HAp'))= 1.91;
inter_a.coupling.scalar(idxof(sys_a,'HA'), idxof(sys_a,'HB')) = 5.05;
inter_a.coupling.scalar(idxof(sys_a,'HA'), idxof(sys_a,'HBp'))= 1.09;
inter_a.coupling.scalar(idxof(sys_a,'HA'), idxof(sys_a,'HC')) =-1.51;
inter_a.coupling.scalar(idxof(sys_a,'HA'), idxof(sys_a,'HCp'))=-1.51;
inter_a.coupling.scalar(idxof(sys_a,'HAp'),idxof(sys_a,'HB')) = 1.09;
inter_a.coupling.scalar(idxof(sys_a,'HAp'),idxof(sys_a,'HBp'))= 5.05;
inter_a.coupling.scalar(idxof(sys_a,'HAp'),idxof(sys_a,'HC')) =-1.51;
inter_a.coupling.scalar(idxof(sys_a,'HAp'),idxof(sys_a,'HCp'))=-1.51;
inter_a.coupling.scalar(idxof(sys_a,'HB'), idxof(sys_a,'HBp'))= 1.98;
inter_a.coupling.scalar(idxof(sys_a,'HB'), idxof(sys_a,'HC')) = 1.33;
inter_a.coupling.scalar(idxof(sys_a,'HB'), idxof(sys_a,'HCp'))= 1.33;
inter_a.coupling.scalar(idxof(sys_a,'HBp'),idxof(sys_a,'HC')) = 1.33;
inter_a.coupling.scalar(idxof(sys_a,'HBp'),idxof(sys_a,'HCp'))= 1.33;
inter_a.coupling.scalar=num2cell(inter_a.coupling.scalar);

% Acrylonitrile (substance B)
props_b=gparse([own_folder 'acrylonitrile.log']);
[sys_b,inter_b]=g2spinach(props_b,{{'H','1H'}},31.8);

% Replace isotropic shifts 
inter_b.zeeman.matrix=shift_iso(inter_b.zeeman.matrix,[1 2 3],...
                                [6.203 6.077 5.703]);

% Create labels 
sys_b.labels={'H7','H8','H9'};

% Replace couplings
inter_b.coupling.scalar=zeros(3,3);
inter_b.coupling.scalar(idxof(sys_b,'H7'),idxof(sys_b,'H8'))= 0.91;
inter_b.coupling.scalar(idxof(sys_b,'H7'),idxof(sys_b,'H9'))=17.92;
inter_b.coupling.scalar(idxof(sys_b,'H8'),idxof(sys_b,'H9'))=11.75;
inter_b.coupling.scalar=num2cell(inter_b.coupling.scalar);

% Endo-norbornene carbonitrile (substance C)
props_c=gparse([own_folder 'norbornene_endo.log']);
[sys_c,inter_c]=g2spinach(props_c,{{'H','1H'}},31.8);

% Replace isotropic shifts 
inter_c.zeeman.matrix=shift_iso(inter_c.zeeman.matrix,[1 2 3 4 5 6 7 8 9],...
                               [6.262 6.129	2.074 1.263	2.781 3.167 ...
                           	    2.959 1.133 1.449]);

% Create labels
sys_c.labels={'H10','H11','H14','H15','H13','H12','H16','H17','H18'};

% Replace couplings (signs to be added)
inter_c.coupling.scalar=zeros(9,9);
inter_c.coupling.scalar(idxof(sys_c,'H12'), idxof(sys_c,'H11'))= 2.8;
inter_c.coupling.scalar(idxof(sys_c,'H12'), idxof(sys_c,'H13'))= 3.8;
inter_c.coupling.scalar(idxof(sys_c,'H13'), idxof(sys_c,'H15'))= 3.8;
inter_c.coupling.scalar(idxof(sys_c,'H14'), idxof(sys_c,'H15'))= 11.8;
inter_c.coupling.scalar(idxof(sys_c,'H14'), idxof(sys_c,'H13'))= 9.42;
inter_c.coupling.scalar(idxof(sys_c,'H14'), idxof(sys_c,'H16'))= 3.4; 
inter_c.coupling.scalar(idxof(sys_c,'H15'), idxof(sys_c,'H16'))= 4.8;
inter_c.coupling.scalar(idxof(sys_c,'H10'), idxof(sys_c,'H11'))= 5.7;
inter_c.coupling.scalar(idxof(sys_c,'H10'), idxof(sys_c,'H16'))= 3.05; 
inter_c.coupling.scalar(idxof(sys_c,'H18'), idxof(sys_c,'H17'))= 9.1;
inter_c.coupling.scalar(idxof(sys_c,'H18'), idxof(sys_c,'H12'))= 1.8;
inter_c.coupling.scalar(idxof(sys_c,'H18'), idxof(sys_c,'H16'))= 2.4;
inter_c.coupling.scalar(idxof(sys_c,'H17'), idxof(sys_c,'H12'))= 1.1;
inter_c.coupling.scalar(idxof(sys_c,'H17'), idxof(sys_c,'H16'))= 1.2;
inter_c.coupling.scalar=num2cell(inter_c.coupling.scalar);

% Exo-norbornene carbonitrile (substance D)
props_d=gparse([own_folder 'norbornene_exo.log']);
[sys_d,inter_d]=g2spinach(props_d,{{'H','1H'}},31.8);

% Replace isotropic shifts 
inter_d.zeeman.matrix=shift_iso(inter_d.zeeman.matrix,[1 2 3 4 5 6 7 8 9],...
                               [6.099 5.9713 1.8972 1.4956 ...
                                3.158 2.9830 1.4958 1.4958 2.113]);

% Create labels
sys_d.labels={'H19','H20','H23','H24','H21','H25','H26','H27','H22'};

% Replace couplings (signs to be added)
inter_d.coupling.scalar=zeros(9,9);
inter_d.coupling.scalar(idxof(sys_d,'H21'), idxof(sys_d,'H25'))= 1.5;   
inter_d.coupling.scalar(idxof(sys_d,'H21'), idxof(sys_d,'H20'))= 3.1;  
inter_d.coupling.scalar(idxof(sys_d,'H21'), idxof(sys_d,'H26'))= 1.05;   
inter_d.coupling.scalar(idxof(sys_d,'H21'), idxof(sys_d,'H27'))= 1.05;
inter_d.coupling.scalar(idxof(sys_d,'H22'), idxof(sys_d,'H23'))= 4.5; 
inter_d.coupling.scalar(idxof(sys_d,'H22'), idxof(sys_d,'H24'))= 8.5; 
inter_d.coupling.scalar(idxof(sys_d,'H22'), idxof(sys_d,'H21'))= 0.85;
inter_d.coupling.scalar(idxof(sys_d,'H23'), idxof(sys_d,'H24'))= 12.4;
inter_d.coupling.scalar(idxof(sys_d,'H23'), idxof(sys_d,'H25'))= 3.5;
inter_d.coupling.scalar(idxof(sys_d,'H24'), idxof(sys_d,'H25'))= 0.8;
inter_d.coupling.scalar(idxof(sys_d,'H25'), idxof(sys_d,'H27'))= 1.5;
inter_d.coupling.scalar(idxof(sys_d,'H25'), idxof(sys_d,'H26'))= 1.5;
inter_d.coupling.scalar(idxof(sys_d,'H19'), idxof(sys_d,'H20'))= 5.7;
inter_d.coupling.scalar(idxof(sys_d,'H19'), idxof(sys_d,'H25'))= 3.1;
inter_d.coupling.scalar=num2cell(inter_d.coupling.scalar);

% Acetonitrile (substance E, solvent)
sys_e.isotopes={'1H','1H','1H'};
sys_e.labels={'H28','H29','H30'};
inter_e.zeeman.matrix={2.0, 2.0, 2.0}*eye(3);
inter_e.coordinates={[]; []; []};
inter_e.coupling.scalar=zeros(3,3);
inter_e.coupling.scalar=num2cell(inter_e.coupling.scalar);

% Merge the spin systems
[sys,inter]=merge_inp({sys_a,   sys_b,   sys_c,   sys_d,   sys_e},...
                      {inter_a, inter_b, inter_c, inter_d, inter_e});

% Chemical parts and unit concentrations
inter.chem.parts={1:6, 7:9, 10:18, 19:27, 28:30};
inter.chem.concs=[1 1 1 1 1];

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='IK-2';
bas.connectivity='scalar_couplings';
bas.space_level=1;

% Relaxation theory parameters
inter.relaxation={'redfield','t1_t2'};
inter.equilibrium='zero';
inter.rlx_keep='secular';
inter.tau_c={ 5e-12 ... % Cyclopentadiene
             20e-12 ... % Acrylonitrile
             50e-12 ... % Exo-norbornene carbonitrile
             50e-12 ... % Endo-norbornene carbonitrile
              1e-12};   % Acetonitrile
inter.r1_rates=cell(30,1); inter.r1_rates(:)={0};
inter.r2_rates=cell(30,1); inter.r2_rates(:)={0};
inter.r1_rates(28:30)={0.5}; % Solvent
inter.r2_rates(28:30)={0.5}; % Solvent

end

% Alles am Weibe ist ein Rätsel, und Alles am Weibe 
% hat Eine Lösung: sie heißt Schwangerschaft.
%
% Friedrich Nietzsche

% #NGRUM