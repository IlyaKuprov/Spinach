% Reproduction of MOLCAS results with the Ligand Field Theory model 
% for a single Dy(III) ion.
%
% Calculation time: seconds
%
% e.suturina@soton.ac.uk
% i.kuprov@soton.ac.uk

function dy_lft_single_1()

% Magnetic field
sys.magnet=0;

% Single Dy ion
sys.isotopes={'E16'};

% Real g-tensor
D=[ 1.325781    1.322640    1.317917];
V=[-0.507708    0.520032    0.686877
   -0.589399    0.371842   -0.717177
   -0.628364   -0.768961    0.117719];
inter.zeeman.matrix={V'*diag(D)*V};

% Rotate the ligand field into the molecular frame
R=[-0.00244652355382  0.10885648156199  0.99405446578366
   -0.90442270078080  0.42385563577788 -0.04864132329304
   -0.42663051090478 -0.89916442681040  0.09741530025544];
[alp,bet,gam]=dcm2euler(R');

% Liza - this needs more decimal places
rkd=[ -0.2821    0.7292   	 0.6234
       0.2751   -0.5610      0.7807
       0.9191    0.3918     -0.0423 ];
[a,b,g]=dcm2euler(rkd');

% Ligand field parameters (MOLCAS)
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
    SBkq{k}=wigner(k,a,b,g)*...
            wigner(k,alp,bet,gam)*stev2sph(k,Bkq{k}); %#ok<AGROW>
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

% Display Spinach answer
[~,D]=eig(geffect(spin_system,[1 2]));
disp('  '); disp('Spinach results:');
disp(['Eigenvalues: ' num2str(real(diag(D)'))]);

% Display MOLCAS answer
disp('  '); disp('MOLCAS results:');
disp('Eigenvalues: 19.2967    0.0529       0.0579');

end

