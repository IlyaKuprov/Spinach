% Accuracy benchmark for the Hamiltonian period propagator
% caclulation using Lie group integrators. Bacause the mo-
% dulation is purely sinusoidal, 2nd and 4th order integra-
% tors show the same apparent accuracy.
%
% Calculation time: seconds
%
% ilya.kuprov@weizmann.ac.il
% u.rasulov@soton.ac.uk

function aht_benchmark()

% System specification
sys.magnet=14.1; sys.isotopes={'14N'};
inter.coupling.matrix{1,1}=eeqq2nqi(1.18e6,0.53,1,[0 0 0]);
inter.zeeman.scalar{1}=32.4;

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Relaxation theory
inter.relaxation={'damp'};
inter.rlx_keep='diagonal';
inter.equilibrium='zero';
inter.damp_rate=300;

% Algorithmic options
sys.disable={'krylov','trajlevel'};

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Magic angle
theta=atan(sqrt(2));

% Spectrum setup
parameters.max_rank=6;
parameters.axis=[sqrt(2/3) 0 sqrt(1/3)];
parameters.rate=-19840;
parameters.grid='single_crystal';
parameters.Lx=cos(theta)*operator(spin_system,'Lz','14N')+...
              sin(theta)*operator(spin_system,'Lx','14N');
parameters.spins={'14N'};
parameters.rf_pwr=2*pi*55e3/sin(theta);
parameters.rf_frq=48e3;

% Pull the Hamiltonian out of the wrapper
stuff=singlerot(spin_system,@impound,parameters,'qnmr'); 
H0=stuff{3}; parameters=stuff{2};

% Get the overtone frequency
ovt_frq=-2*spin('14N')*spin_system.inter.magnet/(2*pi);

% Get the pulse frequency
omega=2*pi*ovt_frq-2*pi*parameters.rf_frq;

% Project pulse operators
Lx=kron(speye(parameters.spc_dim),parameters.Lx);
        
% Get the pulse Hamiltonian
pulseop=parameters.rf_pwr*Lx; Hp=pulseop/2; Hm=pulseop/2;

% Reference calculation using 4th order Lie quadrature
spin_system.tols.prop_chop=eps();
nslices=64; slice_dur=2*pi/(omega*nslices); P_ref=eye(size(H0));
for n=1:nslices
    HL=H0+exp(+1i*omega*slice_dur*(n-1.0))*Hp+exp(-1i*omega*slice_dur*(n-1.0))*Hm;
    HM=H0+exp(+1i*omega*slice_dur*(n-0.5))*Hp+exp(-1i*omega*slice_dur*(n-0.5))*Hm;
    HR=H0+exp(+1i*omega*slice_dur*(n-0.0))*Hp+exp(-1i*omega*slice_dur*(n-0.0))*Hm;
    P_ref=propagator(spin_system,isergen(HL,HM,HR,slice_dur),slice_dur)*P_ref;
end

% Pre-allocate error array
nslices=2:32; error=zeros(numel(nslices),4);

% Loop over point count
for k=1:numel(nslices)

    % Discretise the period of the rotating frame
    slice_dur=2*pi/(omega*nslices(k));

    % Left edge quadrature
    PL=eye(size(H0)); 
    parfor n=1:nslices(k)
        HL=H0+exp(+1i*omega*slice_dur*(n-1.0))*Hp+exp(-1i*omega*slice_dur*(n-1.0))*Hm;
        PL=propagator(spin_system,HL,slice_dur)*PL;
    end
    error(k,1)=norm(P_ref-PL,1)/norm(P_ref,1);

    % Midpoint quadrature
    PM=eye(size(H0)); 
    parfor n=1:nslices(k)
        HM=H0+exp(+1i*omega*slice_dur*(n-0.5))*Hp+exp(-1i*omega*slice_dur*(n-0.5))*Hm;
        PM=propagator(spin_system,HM,slice_dur)*PM;
    end
    error(k,2)=norm(P_ref-PM,1)/norm(P_ref,1);

    % Second order Lie quadrature
    P2=eye(size(H0)); 
    parfor n=1:nslices(k)
        HL=H0+exp(+1i*omega*slice_dur*(n-1.0))*Hp+exp(-1i*omega*slice_dur*(n-1.0))*Hm;
        HR=H0+exp(+1i*omega*slice_dur*(n-0.0))*Hp+exp(-1i*omega*slice_dur*(n-0.0))*Hm;
        P2=propagator(spin_system,isergen(HL,[],HR,slice_dur),slice_dur)*P2;
    end
    error(k,3)=norm(P_ref-P2,1)/norm(P_ref,1);

    % Fourth order Lie quadrature
    P4=eye(size(H0));
    parfor n=1:nslices(k)
        HL=H0+exp(+1i*omega*slice_dur*(n-1.0))*Hp+exp(-1i*omega*slice_dur*(n-1.0))*Hm;
        HM=H0+exp(+1i*omega*slice_dur*(n-0.5))*Hp+exp(-1i*omega*slice_dur*(n-0.5))*Hm;
        HR=H0+exp(+1i*omega*slice_dur*(n-0.0))*Hp+exp(-1i*omega*slice_dur*(n-0.0))*Hm;
        P4=propagator(spin_system,isergen(HL,HM,HR,slice_dur),slice_dur)*P4;
    end
    error(k,4)=norm(P_ref-P4,1)/norm(P_ref,1);

end

% Plotting
figure(); grid on; plot(nslices',error);
set(gca,'YScale','log','XScale','log');
set(gca,'MinorGridLineStyle','-'); 
set(gca,'MinorGridColor',0.9*[1 1 1]);
set(gca,'GridColor',0.9*[1 1 1]);
kxlabel('number of points in the RF period grid');
kylabel('$\|$difference$\|$/$\|$exact$\|$'); 
klegend({'LP','MP','LG-2','LG-4'},...
        'Location','SouthWest');
kgrid; xlim tight; ylim([1e-10 1e-2]);
scale_figure([1.00 0.65]);

end

