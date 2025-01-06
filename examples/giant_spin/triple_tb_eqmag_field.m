% Simulation of the field dependence of the magnetisation of a triple-
% Tb triangular complex - see Figure S27 and S28 in the Supplementary
% Information of the following paper
%
%                https://doi.org/10.1002/chem.201703842
%
% Ligand field parameters and g-tensor for the J=6 ground term were
% computed using the SINGLE_ANISO routine in MOLCAS.
%
% Calculation time: hours
%
% e.suturina@soton.ac.uk
% ilya.kuprov@weizmann.ac.il

function triple_tb_eqmag_field()

% Three J=6 terbium atoms
sys.isotopes={'E13','E13','E13'};

% g-tensor eigenvalues
g1v=[1.497075749 1.495252923 1.481370349];
g2v=[1.496540374 1.494559188 1.482686210];
g3v=[1.497265940 1.494866241 1.481858735];

% g-tensor eigenvalues (rows)
U1=[-0.847221 -0.510119 -0.148306
    -0.119123 -0.089636  0.988825
    -0.517712  0.855420  0.015175];
U2=[ 0.632422 -0.325205  0.703053
    -0.595957  0.375537  0.709794
    -0.494851 -0.867879  0.043689];
U3=[ 0.038754 -0.801228  0.597102
     0.102049 -0.591253 -0.800003
     0.994024  0.091937  0.058852];
 
% g-tensor matrices 
g1=U1*diag(g1v)*U1';
g2=U2*diag(g2v)*U2';
g3=U3*diag(g3v)*U3';
inter.zeeman.matrix={g1,g2,g3};

% Tb atom coordinates
inter.coordinates={[ 2.074806000      1.089908000     -0.432651000]
                   [-2.010910000      1.254299000     -0.348424000]
                   [-0.076869000     -2.348277000     -0.304763000]};
               
% Exchange couplings (Spinach uses NMR convention)
J=icm2hz(0.003);
inter.coupling.scalar={0 J J; 
                       0 0 J;
                       0 0 0};

% Ligand field rotations
R1=[-0.84722146280957 -0.51011886371357 -0.14830555565586
    -0.11912258988979 -0.08963606761341  0.98882515338191
    -0.51771189046879  0.85542043479588  0.01517492012673];
R2=[ 0.63242234003580 -0.32520539385769  0.70305293942172
    -0.59595675112025  0.37553661563385  0.70979419630796
    -0.49485102265009 -0.86787885221622  0.04368939525808];
R3=[ 0.03875378043232 -0.80122835625537  0.59710239124837
     0.10204946466104 -0.59125348643743 -0.80000326345458
     0.99402417036237  0.09193713019178  0.05885161703387];
[alp2,bet2,gam2]=dcm2euler(R2');
[alp1,bet1,gam1]=dcm2euler(R1');
[alp3,bet3,gam3]=dcm2euler(R3');

% Stevens coefficients, Tb 1
Bkq1{2}=[ 0.14292463619745E+00
          0.21658923202324E+00
         -0.33140926967160E+01
          0.60032391844434E-01
         -0.34446225851076E-01];
Bkq1{4}=[ 0.23860615148267E-01
          0.12761921726717E-01
          0.29178580341228E-01
          0.12249490958578E-01
         -0.24899958599662E-02
          0.10883918241352E-01
          0.49198314453769E-01
         -0.52016603964636E-01
         -0.63954337288585E-02];
Bkq1{6}=[ 0.30409866363673E-03
          0.29055196625909E-03
         -0.24291094887073E-03
          0.64765806167618E-04
         -0.29681455189744E-03
         -0.17518004787538E-03
          0.20738256831800E-04
         -0.10671420915421E-03
         -0.33590415849587E-03
          0.49759784841870E-04
          0.22330758826760E-03
         -0.51574662254096E-03
         -0.61867616051731E-05];
     
% Stevens coefficients, Tb 2            
Bkq2{2}=[-0.51308162627882E+00
          0.12361510240242E+00
         -0.29158352524191E+01
         -0.13816191875976E+00
          0.60439671867331E+00];
Bkq2{4}=[-0.11150671922993E-01
         -0.49436450765688E-01
          0.34377773356676E-01
          0.10415820652642E-01
         -0.28467231045610E-02
         -0.23517312106687E-01
         -0.45075159919989E-01
          0.44560867597321E-01
          0.18374163635278E-01];
Bkq2{6}=[ 0.21927537442048E-03
          0.38729833587733E-03
          0.85024371415302E-04
          0.42987699294441E-04
         -0.16454995657660E-03
         -0.71362754056714E-04
          0.20788985051287E-04
          0.15821585623698E-03
          0.39742053551055E-03
         -0.72357337431331E-04
         -0.35652585204448E-03
          0.10321196261628E-03
          0.17006270652148E-03];

% Stevens coefficients, Tb 3
Bkq3{2}=[ 0.60574413436613E+00
         -0.16204396070300E+00
         -0.31067317436176E+01
         -0.27421710020226E+00
          0.65151096146868E+00];
Bkq3{4}=[-0.30653118498582E-02
         -0.12350202833528E+00
         -0.15320664291983E-01
         -0.24104929479275E-01
         -0.36788283598152E-02
         -0.14303634455007E-01
          0.39627670779535E-01
         -0.99068260124337E-03
          0.23074798079080E-01];
Bkq3{6}=[-0.17289907271010E-03 
         -0.15619075817764E-03 
         -0.14202820722595E-03 
          0.38681172246870E-03 
          0.45947107590599E-04 
          0.17542915055988E-03 
          0.23660815730785E-04 
          0.20081385547536E-03 
         -0.33649492646116E-03 
         -0.10978719724285E-03 
         -0.26497047244238E-03 
         -0.93670145366188E-03 
          0.21149957398903E-03];

% Convert to irreducible spherical tensors
for k=2:2:6
    Bkq1{k}=icm2hz(Bkq1{k});
    Bkq2{k}=icm2hz(Bkq2{k});
    Bkq3{k}=icm2hz(Bkq3{k});
    Bkq1{k}=wigner(k,alp1,bet1,gam1)*stev2sph(k,Bkq1{k});
    Bkq2{k}=wigner(k,alp2,bet2,gam2)*stev2sph(k,Bkq2{k});
    Bkq3{k}=wigner(k,alp3,bet3,gam3)*stev2sph(k,Bkq3{k});
end

% Supply to Spinach
inter.giant.coeff={{[0 0 0],Bkq1{2},[0 0 0 0 0 0 0],Bkq1{4},[0 0 0 0 0 0 0 0 0 0 0],Bkq1{6}}...
                   {[0 0 0],Bkq2{2},[0 0 0 0 0 0 0],Bkq2{4},[0 0 0 0 0 0 0 0 0 0 0],Bkq2{6}}...
                   {[0 0 0],Bkq3{2},[0 0 0 0 0 0 0],Bkq3{4},[0 0 0 0 0 0 0 0 0 0 0],Bkq3{6}}};
inter.giant.euler={{[0 0 0],[0 0 0],[0 0 0],[0 0 0],[0 0 0],[0 0 0]}...
                   {[0 0 0],[0 0 0],[0 0 0],[0 0 0],[0 0 0],[0 0 0]}...
                   {[0 0 0],[0 0 0],[0 0 0],[0 0 0],[0 0 0],[0 0 0]}};

% Basis set
bas.formalism='zeeman-hilb';
bas.approximation='none';

% Spherical powder grid
parameters.grid='leb_2ang_rank_11';

% Temperature in Kelvin
inter.temperature=2.0;

% Magnetic field range
B0=[0.01 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 ...
    0.9 1 1.1 1.2 1.3 1.4 1.5 2 3 4 5 6];

% Preallocate the answer
Mz=nan(size(B0)); figure();

% Scan magnetic field
for n=1:numel(B0)
    
    % Set the field
    sys.magnet=B0(n);
    
    % Run Spinach housekeeping
    spin_system=create(sys,inter);
    spin_system=basis(spin_system,bas);
    
    % Compute equilibrium magnetisation
    mag=eqmag(spin_system,parameters);
    
    % Take the Z component
    Mz(n)=mag(3);
    
    % Do the plotting (theory)
    plot(B0,Mz); hold on;
    kxlabel('Magnetic field, Tesla');
    kylabel('Magnetisation, Bohr magneton');
    
    % Do the plotting (experiment)
    load('triple_tb_eqmag.mat','field','magn');
    plot(field,magn,'o'); hold off;
    kgrid; box on; axis tight; drawnow();
    
end

end

