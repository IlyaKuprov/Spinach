% Third-spin-assisted recoupling in the Grenoble PAIN-CP, PAR, and
% PSPAR model systems. PAIN-CP uses the same RF field on carbon and
% nitrogen. For PAR and PSPAR the two carbons share one RF field.
% This reduced matching grid illustrates the signed terminal signal;
% it is not a converged reproduction of the full powder maps. Models
% and terminal-sample conventions: dnp-grenoble source revision
% fa4d4d9ec2f49aed9a01bac7343f9e166121ad34, tsar/pain_grid.in,
% tsar/par_grid.in, and tsar/pspar_grid.in (RF matching maps).
%
% talos@spindynamics.org

function tsar_matching()

% Dimensionless RF axes, deliberately unequal for the axis-order check
hetero_rf=0:1:4;
proton_rf=[0 0.8 1.6 2.4];
sys.parallel={'processes',4};
sys.output='hush';
bas.formalism='zeeman-hilb';
bas.approximation={'none'};
parameters.axis=[sqrt(2/3) 0 sqrt(1/3)];
parameters.max_rank=16;
parameters.grid='leb_3ang_rank_5';
panels={'PAIN-CP: full','PAIN-CP: direct only','PAIN-CP: mediator only',...
        'PAR','PSPAR'};

% Compute the three PAIN-CP mechanisms and two carbon transfer controls
kfigure();
tiledlayout(2,3);
for n=1:numel(panels)

    % Set up the model and preserve the source acquisition intervals
    sys.magnet=2*pi*400e6/spin('1H');
    parameters.rate=20000;
    parameters.samples=0:999;
    parameters.reverse_rf=false;
    if n<=3

        % PAIN-CP interaction tensors in Hz
        sys.isotopes={'13C','15N','1H'};
        couplings=[938.47 -3564.64 12252.69];
        angles=[177.66 49.41 174.90; 171.48 34.70 -177.84; 6.35 15.31 -8.00];
        inter.zeeman.scalar={0 0 0};

        % Isolate the direct pair or the two proton links
        if n==2
            couplings(2:3)=0;
        elseif n==3
            couplings(1)=0;
        end

    else

        % PAR and PSPAR use the same carbon pair and the deck's site shifts
        sys.isotopes={'13C','13C','1H'};
        couplings=[-2139.35 -21521.22 -2955.82];
        angles=[-115.31 50.24 125.94; 113.36 31.72 -124.77; -101.56 72.58 92.23];
        shifts=[-1006 2012 0.430 90 90 180; 1006 2012 0.500 87 64 338;...
                 1000 -3080 0.650 90 90 90];
        inter.zeeman.matrix=cell(1,3);
        for k=1:3
            hz_per_ppm=1e-6*sys.magnet*spin(sys.isotopes{k})/(2*pi);
            euler=deg2rad(shifts(k,4:6));
            inter.zeeman.matrix{k}=-anas2mat(shifts(k,1),1.5*shifts(k,2),...
                                             shifts(k,3),-euler(3),-euler(2),...
                                             -euler(1))/hz_per_ppm;
        end
        inter.zeeman=struct('matrix',{inter.zeeman.matrix});

    end

    % Use the source PSPAR field, MAS rate, and opposite-phase second block
    if n==5
        sys.magnet=2*pi*900e6/spin('1H');
        inter.zeeman.matrix=cellfun(@(x) x*(400/900),inter.zeeman.matrix,...
                                   'UniformOutput',false);
        parameters.rate=60000;
        parameters.reverse_rf=true;
        parameters.samples=2400;
    end

    % Assemble the three dipolar tensors in the active ZYZ convention
    inter.coupling.matrix=cell(3);
    pairs=[1 2; 1 3; 2 3];
    for k=1:3
        Q=euler2dcm(deg2rad(-fliplr(angles(k,:))));
        inter.coupling.matrix{pairs(k,1),pairs(k,2)}=...
            couplings(k)*Q*diag([-1 -1 2])*Q';
    end

    % Spinach housekeeping and spin-resolved preparation and detection
    spin_system=create(sys,inter);
    spin_system=basis(spin_system,bas);
    parameters.spins=unique(sys.isotopes,'stable');
    parameters.offset=zeros(size(parameters.spins));
    parameters.rho0=state(spin_system,'Lx',1);
    parameters.coil=coil_state(spin_system,'Lx',2);
    parameters.heavy_x=operator(spin_system,'Lx',[1 2]);
    parameters.proton_x=operator(spin_system,'Lx',3);
    parameters.heavy_rf=hetero_rf;
    parameters.proton_rf=proton_rf;

    % Average the complex traces before selecting the signed terminal point
    traces=singlerot(spin_system,@matching_traces,parameters,'nmr');
    nexttile;
    imagesc(hetero_rf,proton_rf,real(traces(:,:,end)));
    axis xy;
    colorbar;
    kxlabel('Heavy-spin RF / MAS rate');
    kylabel('Proton RF / MAS rate');
    ktitle(panels{n});

end

end

% Evolve midpoint rotor slices and read the source acquisition ticks
function traces=matching_traces(~,parameters,H,~,~)

% Interpolate the supplied Fourier grid to negative-going slice midpoints
ticks=16;
phase=-2*pi*((1:ticks)-0.5)/ticks;
rank=(numel(H)-1)/2;
shape=size(H{1});
samples=cellfun(@(x) full(x(:)),H,'UniformOutput',false);
coeffs=fftshift(fft([samples{:}],[],2),2)/numel(H);
midpoints=coeffs*exp(1i*(-rank:rank).'*phase);
dwell=1/(ticks*parameters.rate);
traces=complex(zeros(numel(parameters.proton_rf),numel(parameters.heavy_rf),...
                     numel(parameters.samples)));

% Compute one coherent trace for each independently specified RF pair
for n=1:numel(parameters.proton_rf)
    for k=1:numel(parameters.heavy_rf)

        % Form the two-channel RF Hamiltonian in radians per second
        rf=2*pi*parameters.rate*(parameters.heavy_rf(k)*parameters.heavy_x+...
                                 parameters.proton_rf(n)*parameters.proton_x);
        prefix=cell(1,ticks+1);
        prefix{1}=eye(shape);
        for s=1:ticks
            drift=reshape(midpoints(:,s),shape);
            prefix{s+1}=expm(-1i*(drift+rf)*dwell)*prefix{s};
        end

        % PSPAR applies 150 positive periods followed by 150 negative periods
        if parameters.reverse_rf

            % Form the second period with both RF phases inverted
            negative=eye(shape);
            for s=1:ticks
                drift=reshape(midpoints(:,s),shape);
                negative=expm(-1i*(drift-rf)*dwell)*negative;
            end
            U=negative^150*prefix{end}^150;
            traces(n,k,1)=trace(parameters.coil'*U*parameters.rho0*U')/shape(1);

        else

            % Acquire at t=0 and after each source dwell, not at the maximum
            U=eye(shape);
            for s=1:numel(parameters.samples)
                tick=parameters.samples(s);
                if tick>0&&mod(tick,ticks)==0
                    U=prefix{end}*U;
                end
                P=prefix{mod(tick,ticks)+1}*U;
                traces(n,k,s)=trace(parameters.coil'*P*parameters.rho0*P')/shape(1);
            end

        end

    end
end

end


