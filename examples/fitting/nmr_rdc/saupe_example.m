% Extracting Saupe order matrix from NH RDC data. Experimental 
% measurements kindly provided by Andras Boeszoermenyi, Thibault
% Viennet, and Hari Arthanari.
%
% ilya.kuprov@weizmann.ac.il

function saupe_example()

% Read the PDB file
[pdb_aa,~,pdb_id,coords]=read_pdb_pro('protein.pdb',1);

% Read RDC data
load('rdc_data.mat','aa_num','rdc_hz',...
     'isotope_a','isotope_b','atom_a','atom_b');

% Make isotope table
isotopes=[isotope_a isotope_b];

% Match up RDCs with coordinates
for n=1:numel(rdc_hz)
    
    % Locate both atoms
    index_a=(pdb_aa==aa_num(n))&strcmp(pdb_id,atom_a{n});
    index_b=(pdb_aa==aa_num(n))&strcmp(pdb_id,atom_b{n});

    % Extract coordinates
    xyz{n,1}=coords{index_a}; xyz{n,2}=coords{index_b}; %#ok<AGROW> 

end

% Call RDC fitter
S=rdc_fit(isotopes,xyz,rdc_hz); 
disp('Saupe order matrix:'); disp(S);

% Back-calculate RDCs
for n=1:numel(rdc_hz)
    rdc_theo(n)=xyz2rdc(isotope_a{n},isotope_b{n},...
                        xyz{n,1},xyz{n,2},{S,'saupe'}); %#ok<AGROW> 
end

% Do the plotting
figure(); plot(rdc_hz,rdc_theo','r.'); 
kxlabel('Experimental RDC, Hz'); 
kylabel('Theoretical RDC, Hz');
xlim([-30 30]); ylim([-30 30]); 
axis square; kgrid; hold on;
plot([-30 30],[-30 30],'b-');

end

