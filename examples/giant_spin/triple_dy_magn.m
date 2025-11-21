% Simulation of a finite-speed magnetic field sweep experiment for a
% single crystal of a triple-Dy triangular complex in a micro-SQUID, 
% see Figure S24 in the Supplementary Information of
%
%                https://doi.org/10.1002/chem.201703842
%
% Ligand field parameters and g-tensor for the J=15/2 ground term were
% computed using the SINGLE_ANISO routine in MOLCAS.
%
% Calculation time: hours
%
% e.suturina@soton.ac.uk
% ilya.kuprov@weizmann.ac.il

function triple_dy_magn()

% Three J=15/2 dysprosium atoms
sys.isotopes={'E16','E16','E16'};

% g-tensor eigenvalues
g_eigs=[1.325781502 1.322640525 1.317917615];

% g-tensor eigenvectors
U3=[ -0.507708  0.520032  0.686877
     -0.589399  0.371842 -0.717177
     -0.628364 -0.768961  0.117719];

% g-tensor matrix
g3=U3'*diag(g_eigs)*U3;

% Triangle arrangement
R120=euler2dcm(0,0,-2*pi/3); 
g2=R120*g3*R120'; g1=R120*g2*R120';

% Data supplied to Spinach
inter.zeeman.matrix={g1,g2,g3};

% Coordinates
inter.coordinates={[ 0.011931000    7.936473000  12.476785000]
                   [ 2.024747000   11.464103000  12.476785000]
                   [-2.036678000   11.443438000  12.476785000]};

% Exchange couplings (Spinach uses NMR convention)
J=icm2hz(0.0063);
inter.coupling.scalar={0 J J; 
                       0 0 J;
                       0 0 0};
                   
% Rotate the ligand field into the molecular frame
R=[-0.00244652355382  0.10885648156199  0.99405446578366
   -0.90442270078080  0.42385563577788 -0.04864132329304
   -0.42663051090478 -0.89916442681040  0.09741530025544];
[alp,bet,gam]=dcm2euler(R');
                   
% Stevens coefficients
Bkq{2}=[ 0.10326442519536E+01
         0.87708188068863E+00
        -0.18677405749006E+01
         0.91671003105747E-01
         0.21762357826531E+00];
Bkq{4}=[ 0.11295055284297E-01
        -0.10066233934224E-01
        -0.13233111716614E-01
        -0.71197744039963E-02
         0.96641829711565E-03
        -0.35617075647864E-02
         0.17555195686494E-01
        -0.14706928649860E-01
         0.83452445804378E-02];
Bkq{6}=[-0.21293172047501E-03
         0.40418140640703E-03
        -0.31028761183788E-04
        -0.53546889849400E-04
         0.10980961975085E-03
         0.16627815256680E-03
        -0.64316087884035E-05
         0.10520043927213E-03
        -0.16007101746115E-03
        -0.24045075061691E-04
        -0.24385297781429E-03
         0.10725179140414E-03
        -0.15569612484251E-03];

% Convert to irreducible spherical tensors
for k=2:2:6
    Bkq{k}=icm2hz(Bkq{k});
    Bkq3{k}=wigner(k,alp,bet,gam)*stev2sph(k,Bkq{k}); %#ok<AGROW>
end

% Supply to Spinach
inter.giant.coeff={{[0 0 0],Bkq3{2},[0 0 0 0 0 0 0],Bkq3{4},[0 0 0 0 0 0 0 0 0 0 0],Bkq3{6}}...
                   {[0 0 0],Bkq3{2},[0 0 0 0 0 0 0],Bkq3{4},[0 0 0 0 0 0 0 0 0 0 0],Bkq3{6}}...
                   {[0 0 0],Bkq3{2},[0 0 0 0 0 0 0],Bkq3{4},[0 0 0 0 0 0 0 0 0 0 0],Bkq3{6}}};
inter.giant.euler={{[+2*pi/3 0 0],[+2*pi/3 0 0],[+2*pi/3 0 0],[+2*pi/3 0 0],[+2*pi/3 0 0],[+2*pi/3 0 0]}...
                   {[-2*pi/3 0 0],[-2*pi/3 0 0],[-2*pi/3 0 0],[-2*pi/3 0 0],[-2*pi/3 0 0],[-2*pi/3 0 0]}...
                   {[0 0 0],[0 0 0],[0 0 0],[0 0 0],[0 0 0],[0 0 0]}};

% Basis set
bas.formalism='zeeman-hilb';
bas.approximation='none';

% Temperature in Kelvin
inter.temperature=0.03;

% This must be set to 1 Tesla
sys.magnet=1.0;

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Experiment parameters
parameters.fields=[0 1];            % Tesla
parameters.npoints=5000;
parameters.sweep_time=1e-5;         % seconds
parameters.orientation=[0 pi/2 0];
parameters.nstates=64;

% Run the experiment
[fields,z_magn]=fieldscan_magn(spin_system,parameters);

% Plot the results
kfigure(); plot(fields,z_magn); 
kxlabel('Magnetic field, Tesla');
kylabel('Magnetisation'); kgrid;

end

