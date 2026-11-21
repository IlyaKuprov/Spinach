% Hahn echo of the negatively charged silicon vacancy in 3C-SiC.
% Fazio et al., doi:10.1103/PhysRevA.109.022603, Figure 12: CCE2 on the
% m_s=+3/2 to +1/2 transition at 200 G. Independent natural-abundance
% baths use ideal zinc-blende geometry, the printed Table II Azz values
% for the first two shells, and point dipoles elsewhere. Unknown transverse terms
% in the hyperfine z row and the original isotope realisations are not recovered.
% The electron Zeeman and axial D=0.1 MHz phases cancel in the ideal echo.
% This reconstruction does not reproduce the published trace quantitatively.
% Calculation time: 291 seconds with 64 process workers on Elminster.
%
%             [answer,diagnostics]=gcce_highspin_siv()
%
% answer is the mean absolute normalised echo, a column over total time
% in diagnostics.times (seconds). Diagnostics also include individual
% bath curves, retained singleton/pair counts, propagated absolute bounds,
% and CCE validity reasons per realisation.
%
% ilya.kuprov@weizmann.ac.il

function [answer,diagnostics]=gcce_highspin_siv()

% Field along [001] and spin-3/2 centre
sys.magnet=0.02;
sys.isotopes={'E4'};
sys.output='hush';
sys.parallel={'processes',8};

% Zinc-blende site catalogue centred at a silicon vacancy
lattice=4.3596;
bath_radius=50;
pair_radius=8;
extent=ceil(bath_radius/lattice);
[x,y,z]=ndgrid(-extent:extent);
fcc=[0 0 0;0 0.5 0.5;0.5 0 0.5;0.5 0.5 0];
si=lattice*reshape(reshape([x(:) y(:) z(:)],[],1,3)+reshape(fcc,1,[],3),[],3);
carbon=si+lattice/4;
xyz=[si;carbon];
species=[ones(size(si,1),1);2*ones(size(carbon,1),1)];
distance=vecnorm(xyz,2,2);
keep=(distance>0)&(distance<=bath_radius);
xyz=xyz(keep,:);
distance=distance(keep);
species=species(keep);

% Natural abundance and independent, reproducible site occupancy
labels={'29Si','13C'};
abundance=[0.047;0.011];
stream=RandStream('mt19937ar','Seed',2024);
nbaths=50;
occupied=rand(stream,numel(species),nbaths)<abundance(species);

% The ideal nearest silicon and carbon distances identify Table II sites
shell_dist=lattice*[sqrt(0.5);sqrt(3)/4];
near_shell=abs(distance-shell_dist(species))<eps(lattice)*4;
table_azz=1e6*[6.8;79.7];

% Bath gyromagnetic ratios and the dipolar prefactor for distances in Angstrom
spin_system=create(sys,struct());
gamma=[spin('29Si');spin('13C')];
dipolar_hz=spin_system.tols.hbar*spin_system.tols.mu0*1e30/(8*pi^2);

% Pure-dephasing hyperfine columns in the Spinach pair-input convention
ort=xyz./distance;
hyperfine=(dipolar_hz*spin('E4')*gamma(species)./distance.^3).*...
          ([0 0 1]-3*ort(:,3).*ort);
hyperfine(near_shell,:)=table_azz(species(near_shell)).*[0 0 1];

% Projected CCE2 with ideal midpoint transition swap and a mixed nuclear bath
parameters.observed=1;
parameters.transition=struct('spins',1,'ms',[3/2 1/2]);
parameters.times=linspace(0,2e-3,801);
parameters.cce.approximation='projected';
parameters.cce.order=2;
parameters.cce.graph=struct('mode','coupling','threshold',0);
parameters.cce.bath_state.mode='mixed_product';
parameters.cce.basis.formalism='zeeman-hilb';
parameters.cce.numerics=struct('atol',0.02,'rtol',0,...
                               'max_memory_bytes',2^33);

% Partition per-bath inputs before distributing the CCE calculations
bath_positions=cell(1,nbaths);
bath_isotopes=cell(1,nbaths);
bath_gammas=cell(1,nbaths);
bath_hfc=cell(1,nbaths);
for n=1:nbaths

    % Materialise the occupied isotope sites and their local couplings
    sites=occupied(:,n);
    bath_positions{n}=xyz(sites,:);
    bath_isotopes{n}=labels(species(sites));
    bath_gammas{n}=gamma(species(sites));
    bath_hfc{n}=hyperfine(sites,:);

end

% Calculate the echo separately for each isotope configuration
curves=zeros(numel(parameters.times),nbaths);
counts=zeros(nbaths,2);
bounds=zeros(numel(parameters.times),nbaths);
reasons=cell(1,nbaths);
parfor n=1:nbaths

    % Select the occupied isotope sites
    bath_xyz=bath_positions{n};
    bath_gamma=bath_gammas{n};
    bath_hyperfine=bath_hfc{n};
    nspins=size(bath_xyz,1)+1;
    bath_sys=sys;
    bath_sys.isotopes=[sys.isotopes bath_isotopes{n}];

    % Work in the refocused electron frame and retain the nuclear Zeeman terms
    inter=struct();
    inter.coupling.matrix=cell(nspins);
    inter.zeeman.matrix=repmat({zeros(3)},1,nspins);
    for k=1:nspins-1

        % The published Hamiltonian retains only the electron z operator
        inter.coupling.matrix{1,k+1}=[zeros(3,2) bath_hyperfine(k,:).'];

        % Add nuclear dipoles only for neighbours within the pair radius
        delta=bath_xyz(k+1:end,:)-bath_xyz(k,:);
        radii=vecnorm(delta,2,2);
        for q=find(radii<=pair_radius).'

            % Convert the pair dipolar tensor to Hz without double counting
            unit=delta(q,:)/radii(q);
            scale=dipolar_hz*bath_gamma(k)*bath_gamma(k+q)/radii(q)^3;
            inter.coupling.matrix{k+1,k+q+1}=scale*(eye(3)-3*(unit.'*unit));

        end

    end

    % Run the material model through the native Hahn-echo context
    bath_system=create(bath_sys,inter);
    [echo,detail]=clexpn(bath_system,@hahn_echo,parameters,'labframe');
    curves(:,n)=abs(echo);
    bounds(:,n)=detail.reduction.abs_error_bound;
    reasons{n}=detail.reduction.samples{1}.reason;
    counts(n,:)=[nspins-1 numel(detail.clusters)-nspins+1];

end

% Plot the isotope-ensemble mean and return its component curves
answer=mean(curves,2);
diagnostics=struct('times',parameters.times,'curves',curves,'counts',counts,...
                   'abs_error_bounds',bounds,'reasons',{reasons});
kfigure();
plot(1e3*parameters.times,answer);
kxlabel('total echo time / ms');
kylabel('coherence magnitude');
kgrid();

end


