% Direct and relayed longitudinal transfer in a model silicon
% network, using the Grenoble SR26 and S3 trains. The direct-only
% and relay-only models share the same shifts and crystal orientation.
% This single-crystal example uses DQ selection between identical arms.
% Pulse schedules: dnp-grenoble revision fa4d4d9ec2f49aed9a01bac7343f9e166121ad34,
% homonuclear_recoupling/SR26_CC.in (fumarate, cf. Teymoori et al.,
% J. Magn. Reson. 261, 205-220, Fig. 3c, PMID 26515279) and
% homonuclear_recoupling/S3SR26_SiSi_varydist.in (silicon DQ buildup).
% This model adds a relay spin and does not replicate either powder result.
%
% talos@spindynamics.org

function sr26_relay()

% Model shifts and dipolar geometry
sys.magnet=2*pi*400e6/spin('1H');
sys.isotopes={'29Si','29Si','29Si'};
sys.parallel={'processes',2};
inter.zeeman.scalar={-60 60 0};
bas.formalism='sphten-liouv';
bas.approximation={'none'};
positions=[0 0 0; 4.2 0 0; 2.0 2.1 0.4];

% Use the source time convention of two recoupling arms per increment
parameters.rate=10000;
parameters.axis=[sqrt(2/3) 0 sqrt(1/3)];
parameters.grid='single_crystal';
parameters.max_rank=6;
parameters.spins={'29Si'};
parameters.counts=0:6;
time=32e3*parameters.counts/parameters.rate;

% Generate and verify all SR26 events before compiling any products
[phase,rf,width]=rsequence(4,26,11,1,parameters.rate,...
                           '90270_pulse','homo_double_quantum_nupicycle');
angle=pi*11/26;
unit=[angle angle+pi -angle -angle+pi].';
block=repmat(unit,13,1);
literal=[block; -block; pi-block; pi+block];
assert(max(abs(exp(1i*phase)-exp(1i*literal)))<32*eps);
assert(max(abs(rf*width-[pi/2 3*pi/2]))<32*eps);
block_time=16/parameters.rate;
assert(abs(numel(phase)*sum(width)/2-block_time)<32*eps(block_time));

% Form the S3 source's four 360/270/90 composite elements
s3_phase=deg2rad([270 450 270 90 270 90 450 270 450 270 90 270]).';
s3_turns=repmat([1 3/4 1/4],1,4);
assert(max(abs(exp(1i*s3_phase).'-[-1i 1i -1i 1i -1i 1i 1i -1i 1i -1i 1i -1i]))<32*eps);
assert(sum(s3_turns)==8);

% Compare both pulse trains with only the direct link or both relay links
curves=zeros(4,numel(time));
for n=1:2

    % Retain only the specified physical dipolar links
    inter.coupling.matrix=cell(3);
    pairs=[1 2; 1 3; 2 3];
    for k=1:3
        if (n==1&&k==1)||(n==2&&k>1)
            [~,~,~,~,Q]=xyz2dd(positions(pairs(k,1),:),positions(pairs(k,2),:),...
                                '29Si','29Si');
            inter.coupling.matrix{pairs(k,1),pairs(k,2)}=Q/(2*pi);
        end
    end

    % Spinach housekeeping and spin-resolved detection
    spin_system=create(sys,inter);
    spin_system=basis(spin_system,bas);
    parameters.rho0=state(spin_system,'Lz',1);
    parameters.coil=-coil_state(spin_system,'Lz',2);

    % Propagate using the verified SR26 and S3 schedules
    parameters.isotope='29Si';
    parameters.sr26=true;
    parameters.phases=phase;
    parameters.rf=rf;
    parameters.widths=width;
    curves(n,:)=singlerot(spin_system,@relay_echo,parameters,'nmr');
    parameters.sr26=false;
    parameters.phases=s3_phase;
    parameters.rf=pi*parameters.rate;
    parameters.widths=2*s3_turns/parameters.rate;
    curves(n+2,:)=singlerot(spin_system,@relay_echo,parameters,'nmr');

end

% Plot the signed silicon transfer
kfigure(); plot(time,real(curves).','o-');
kxlabel('Total recoupling duration, ms');
kylabel('Longitudinal transfer'); kgrid;
klegend({'SR26 pair','SR26 relay','S3 pair','S3 relay'});

end


