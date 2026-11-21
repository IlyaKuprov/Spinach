% Test: kernel/tedor | Finite-pulse TEDOR
% Tests finite-pulse TEDOR with independent Hilbert matrix products.
% Syntax:
%
%                       result=test_tedor()
%
% Outputs:
%
%    result  - complex FID, selection, and null checks
%
% talos@spindynamics.org

function result=test_tedor()

% Register the independent event and pathway checks
result=new_test_result('kernel/tedor','Finite-pulse TEDOR',...
                       'Direct and ZF schedules must retain signed complex signals.');

% Define a detuned heteronuclear pair with a known secular drift
sys.magnet=9.4;
sys.isotopes={'15N','13C'};
inter.coupling.scalar={[] 1350; [] []};
bas.formalism='sphten-liouv';
bas.approximation={'none'};
spin_system=test_spin_system(sys,inter,bas);
parameters.spins=sys.isotopes;
parameters.rate=8000;
parameters.axis=[sqrt(2/3) 0 sqrt(1/3)];
parameters.max_rank=0;
parameters.grid='single_crystal';
parameters.serial=true;
parameters.offset=[355 -215];
parameters.ncycles=[1 0 3];
parameters.npoints=4;
parameters.rf=25e3;

% Assemble the four-level reference independently of Spinach
Ix=kron([0 1; 1 0]/2,eye(2));
Iz=kron(diag([1 -1])/2,eye(2));
Sx=kron(eye(2),[0 1; 1 0]/2);
Sy=kron(eye(2),[0 -1i; 1i 0]/2);
Sz=kron(eye(2),diag([1 -1])/2);
H=2*pi*(1350*Iz*Sz+355*Iz-215*Sz);
drive=2*pi*parameters.rf;
width=0.5/parameters.rf;
period=1/parameters.rate;
delay=@(t) expm(-1i*H*t);
gap=period/2-width;
Ipulse=expm(-1i*(H+drive*Ix)*width);
Xpulse=expm(-1i*(H+drive*Sx)*width);
Ypulse=expm(-1i*(H+drive*Sy)*width);
rot_i=expm(-1i*pi*Ix/2);
rot_s=expm(-1i*pi*Sx/2);

% Calculate the first direct point and stored-propagator acquisition
mix=delay(gap/2)*Xpulse*delay(gap)*Xpulse*delay(gap/2);
receive=delay(gap/2)*Ipulse*delay(gap)*Ipulse*delay(gap/2);
convert=expm(-1i*(H+drive*(Ix+Sx))*width/2);
receiver=Sx+1i*Sy;
ref=zeros(parameters.npoints,numel(parameters.ncycles));
for n=1:numel(parameters.ncycles)

    % Read the signal before each acquisition block
    for k=1:parameters.npoints
        U=receive^(k-1)*convert*mix^parameters.ncycles(n);
        ref(k,n)=trace(receiver'*U*Ix*U')/4;
    end

end
parameters.sequence='direct';
parameters.rho0=state(spin_system,'Lx','15N');
parameters.coil=coil_state(spin_system,'L+','13C');
signal=singlerot(spin_system,@tedor,parameters,'nmr');
result=test_close(result,'direct complex FID',signal,ref,1e-10,0,...
                  'Finite 90-degree conversion and readout use the stored sender inversion train.');


% Construct both ZF paths on the same complete basis
first=Ypulse*delay(gap)*Xpulse*delay(gap);
second=delay(gap)*Ypulse*delay(gap)*Xpulse;
echo=delay(period-width/2)*Ipulse*delay(period-width/2);
mask=zeros(4);
mask([3 8 9 14])=1;
ref=zeros(2,numel(parameters.ncycles));
rejected=zeros(size(parameters.ncycles));
for n=1:numel(parameters.ncycles)

    % Apply individual quarter turns, rather than collapsed production rotations
    count=parameters.ncycles(n);
    U=second^count*echo*first^count;
    prepared=U*Ix*U';
    minus=rot_i*rot_s'*rot_i*rot_s;
    plus=rot_i*rot_s*rot_i*rot_s;
    rho_minus=minus*prepared*minus';
    rho_plus=plus*prepared*plus';
    filtered=mask.*rho_minus;
    ref(:,n)=[trace(Ix*U*(rho_plus-rho_minus)*U');...
              trace(Ix*U*filtered*U')]/4;
    rejected(n)=norm(rho_minus-filtered,'fro');

end
parameters.coil=coil_state(spin_system,'Lx','15N');
parameters.sequence='zf-cycle';
cycled=singlerot(spin_system,@tedor,parameters,'nmr');
parameters.sequence='zf-filter';
filtered=singlerot(spin_system,@tedor,parameters,'nmr');
result=test_close(result,'signed cycle and element filter',[cycled; filtered],ref,1e-10,0,...
                  'Signed receiver cycling and element selection use different operations.');
result=test_true(result,'nontrivial rejected elements',max(rejected)>1e-3,...
                 'The finite-pulse fixture populates coherence outside the mask.');

% Resolve one-cycle dipolar transfer in a tilted single crystallite
inter.coupling.scalar=cell(2);
inter.coupling.eigs={[] [1000 1000 -2000]; [] []};
inter.coupling.euler={[] [0.23 0.68 0]; [] []};
spin_system=test_spin_system(sys,inter,bas);
parameters.max_rank=8;
parameters.sequence='direct';
parameters.coil=coil_state(spin_system,'L+','13C');
dipolar=singlerot(spin_system,@tedor,parameters,'nmr');
result=test_true(result,'resolved direct transfer',norm(dipolar(:,1))>1e-3,...
                 'Rotor-synchronised inversion recouples a tilted dipole.');

% Turn off the dipole without removing RF or detuning
inter.coupling.eigs{1,2}=[];
inter.coupling.euler{1,2}=[];
spin_system=test_spin_system(sys,inter,bas);
parameters.sequence='direct';
parameters.coil=coil_state(spin_system,'L+','13C');
uncoupled=singlerot(spin_system,@tedor,parameters,'nmr');
result=test_close(result,'uncoupled direct transfer',uncoupled,zeros(size(ref,2),parameters.npoints).',1e-10,0,...
                  'Independent spins cannot transfer magnetisation.');

end


