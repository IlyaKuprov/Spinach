% A demonstration that most lanthanide complexes are in the 
% ZFS limit for the purposes of relaxation theory. One of the
% figures from our forthcoming papers on the subject.
%
% Calculation time: hours.
%
% e.suturina@soton.ac.uk
% i.kuprov@soton.ac.uk

function dy_lft_single_2()

% Magnetic field
sys.magnet=1.0;

% Single Dy ion
sys.isotopes={'E16'};

% Real g-tensor
D=[ 1.322766699  1.324261429  1.328750739];
V=[ 0.566995 -0.819250  0.085715
    0.762814  0.561489  0.320696
   -0.310858 -0.116448  0.943296];
inter.zeeman.matrix={V'*diag(D)*V};

% Rotate the ligand field into the molecular frame
R=[ 0.56699461681939 -0.81924972266747  0.08571462189795
    0.76281364538368  0.56148878891791  0.32069562257063
   -0.31085759909370 -0.11644840844243  0.94329598814843];
[alp,bet,gam]=dcm2euler(R');

% Ligand field parameters (MOLCAS)
Bkq{2}=[ 0.14365204369342E-01
        -0.61531135748287E-01
         0.23619823798400E+01
        -0.15447162033345E+00
        -0.10971445533015E+01];
Bkq{4}=[-0.32862425456238E-01
         0.24378068363721E-01
        -0.57191655094048E-02
        -0.27467726145967E-01
         0.31055525736595E-02
         0.72000950135235E-02
        -0.13083689485342E-01
        -0.32573249739364E-01
        -0.16118273163304E-01];
Bkq{6}=[-0.61447247154955E-04
         0.62131530851511E-03
         0.14999995093883E-03
         0.10986805603157E-03
        -0.16773172829085E-03
        -0.17727850671693E-03
         0.62056810317789E-05
        -0.25339852925514E-04
        -0.11398439843870E-03
         0.38226216128378E-03
        -0.28082485105398E-05
         0.65043175487855E-03
         0.89676120791658E-04];
     
% Convert to irreducible spherical tensors
for k=[2 4 6]
    Bkq{k}=icm2hz(Bkq{k});
    SBkq{k}=wigner(k,alp,bet,gam)*stev2sph(k,Bkq{k}); %#ok<AGROW>
end

% Supply to Spinach
inter.giant.coeff={{[0 0 0],SBkq{2},...
                    [0 0 0 0 0 0 0],SBkq{4},...
                    [0 0 0 0 0 0 0 0 0 0 0],SBkq{6}}};
inter.giant.euler={{[0 0 0],[0 0 0],[0 0 0],...
                    [0 0 0],[0 0 0],[0 0 0]}};

% Formalism specification
bas.formalism='zeeman-hilb';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Experiment parameters
parameters.fields=[0 500];
parameters.npoints=1000;
parameters.orientation=[0 0 0];
parameters.nstates=16;

% Run the simulation
fieldscan_enlev(spin_system,parameters);

end

