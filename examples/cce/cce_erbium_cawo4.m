% Reconstructs the lattice and amorphous Er:CaWO4 Hahn echoes in Fig. S12
% of Le Dantec et al., doi:10.1126/sciadv.abj9786. The curves are
% ensemble means over 50 independently sampled baths of each type.
% Numerical inputs are from the authors' Fortran and supplemental theory
% in doi:10.7910/DVN/9VSAUT. Original bath realisations and seeds were not
% archived; the continuum is sampled uniformly with no hard core.
%
% ilya.kuprov@weizmann.ac.il

function cce_erbium_cawo4()

% Experimental geometry and total Hahn echo times
cell_size=1e-10*[5.2429 5.2429 11.3737];
radius=11e-9;
abundance=0.145;
parameters.observed=1;
parameters.transition=struct('spins',1,'ms',[1/2 -1/2]);
parameters.orientation=deg2rad([0 90 133.5]);
parameters.times=(0:60)*1e-3;
parameters.cce.approximation='projected';
parameters.cce.order=2;
parameters.cce.numerics=struct('atol',1e-4,'max_memory_bytes',2^27);
n_configs=50;

% Four W sites relative to a Ca site in the authors' conventional cell
sites=[0.0 0.5 -0.25;...
       0.5 0.5  0.00;...
       0.5 0.0  0.25;...
       0.0 0.0 -0.50];
extent=ceil(radius./cell_size);
[ax,ay,az]=ndgrid(-extent(1):extent(1),...
                 -extent(2):extent(2),...
                 -extent(3):extent(3));
cells=[ax(:) ay(:) az(:)];
lattice=(repmat(cells,4,1)+repelem(sites,size(cells,1),1)).*cell_size;
lattice=lattice(vecnorm(lattice,2,2)<=radius,:);
n_amorphous=10346;

% Express the paper's dipolar prefactors in Hz Angstrom cubed
mu0=1.25663706127e-6;
hbar=6.62607015e-34/(2*pi);
mu_bohr=9.2740100657e-24;
gamma_w=spin('183W');
g_perp=8.38;
convert=1e30/(2*pi);
inter.dipolar_geometry.hyperfine=convert*mu0/(4*pi)*g_perp*mu_bohr*gamma_w;
inter.dipolar_geometry.nuclear=convert*mu0/(4*pi)*hbar*gamma_w^2;
sys.magnet=0.067;
sys.parallel={'processes',8};
sys.output='hush';

% Prepare independent isotope occupancies and continuum configurations
occupied=cell(n_configs,1);
amorphous=cell(n_configs,1);
for n=1:n_configs

    % Select isotope sites by Bernoulli occupancy and sample uniform continuum positions
    stream=RandStream('mrg32k3a','Seed',1);
    stream.Substream=n;
    occupied{n}=lattice(rand(stream,size(lattice,1),1)<abundance,:);
    directions=randn(stream,n_amorphous,3);
    amorphous{n}=radius*rand(stream,n_amorphous,1).^(1/3).*directions./vecnorm(directions,2,2);

end

% Set up the first physical bath on the head node so create configures the pool
first_bath=bath_system(sys,inter,occupied{1});

% Evaluate each independently created bath on a worker
realisations=zeros(numel(parameters.times),n_configs,2);
parfor n=1:n_configs

    % Reuse the first created system instead of creating it again
    if n==1
        lattice_bath=first_bath;
    else
        lattice_bath=bath_system(sys,inter,occupied{n});
    end
    continuum_bath=bath_system(sys,inter,amorphous{n});
    lattice_echo=clexpn(lattice_bath,@hahn_echo,parameters,'labframe');
    random_echo=clexpn(continuum_bath,@hahn_echo,parameters,'labframe');
    realisations(:,n,:)=reshape([lattice_echo random_echo],numel(parameters.times),1,2);

end
echo=reshape(mean(realisations,2),numel(parameters.times),2);

% Plot the normalised amplitude, not its square
kfigure();
plot(1e3*parameters.times,echo,'LineWidth',1.5);
kxlabel('total echo time, $2\tau$ / ms');
kylabel('normalised Hahn echo amplitude');
klegend({'lattice bath','amorphous bath'},'Location','northeast');
kgrid();

end

% Create one compact physical bath with the central electron first
function spin_system=bath_system(sys,inter,coords)

% Assign physical IDs with the central electron first
sys.isotopes=[{'E'} repmat({'183W'},1,size(coords,1))];
inter.coordinates=num2cell([zeros(1,3); 1e10*coords],2);
spin_system=create(sys,inter);

end


