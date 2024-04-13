% Nuclear relaxation rates using the adiabatic elimination method
% for a rapidly relaxing Dy(III) ion with a user-specified ZFS.
%
% Calculation time: minutes
%
% e.suturina@soton.ac.uk
% i.kuprov@soton.ac.uk

function nuclear_relaxation_1()

% Magnetic field
sys.magnet=14.1;

% Dy(III) ion and a proton
sys.isotopes={'E16','1H'};

% Electron g-tensor
D=[ 1.325781    1.322640    1.317917];
V=[-0.507708    0.520032    0.686877
   -0.589399    0.371842   -0.717177
   -0.628364   -0.768961    0.117719];
inter.zeeman.matrix{1}=V'*diag(D)*V;

% Nuclear shift tensor
inter.zeeman.matrix{2}=zeros(3,3);

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
                    [0 0 0 0 0 0 0 0 0 0 0],SBkq{6}},{}};
inter.giant.euler={{[0 0 0],[0 0 0],[0 0 0],...
                    [0 0 0],[0 0 0],[0 0 0]},{}};

% Cartesian coodinates
inter.coordinates={[0.00 0.00 0.00];
                   [0.00 5.00 7.00]};
               
% Separate spin relaxation times
inter.relaxation={'t1_t2'};
inter.r1_rates={1/50e-15 0}; % T1e = 50 fs
inter.r2_rates={1/50e-15 0}; % T2e = 50 fs
inter.rlx_keep='labframe';
inter.equilibrium='zero';

% Formalism specification
bas.formalism='sphten-liouv';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Slow and fast subspace indices
slow_idx=1:4;    % First four states are pure nuclear
fast_idx=5:1024; % Other states involve the electron

% Hamiltonian
[I,Q]=hamiltonian(assume(spin_system,'labframe'));

% Non-interacting relaxation
R_ni=relaxation(spin_system);

% Get spherical averaging grid as a structure
sph_grid=load([spin_system.sys.root_dir '/kernel/grids/' ...
               'leb_2ang_rank_11.mat'],'alphas','betas',...
               'gammas','weights');
alphas=sph_grid.alphas; betas=sph_grid.betas;
gammas=sph_grid.gammas; weights=sph_grid.weights;

% Orientation average (3 angles)
R=zeros(4,4,'like',1i);
parfor n=1:numel(weights)
    
    % Liouvillian with non-interacting relaxation
    L=I+orientation(Q,[alphas(n) betas(n) gammas(n)])+1i*R_ni;

    % Adiabatic elimination at the current orientation
    [~,R_curr]=adelim(spin_system,L,fast_idx,slow_idx);
    
    % Spherical average
    R=R+weights(n)*R_curr;
    
end

% Get the observables in the nuclear subspace
Nz=state(spin_system,{'Lz'},{2}); Nz=Nz(slow_idx);
Np=state(spin_system,{'L+'},{2}); Np=Np(slow_idx);
Nz=Nz/norm(Nz,2); Np=Np/norm(Np,2);

% Report relaxation rates
disp(['1H R1:  ' num2str(real(-Nz'*R*Nz)) ' Hz']);
disp(['1H R2:  ' num2str(real(-Np'*R*Np)) ' Hz']);
disp(['1H DFS: ' num2str(imag(-Np'*R*Np)) ' Hz']);

end

