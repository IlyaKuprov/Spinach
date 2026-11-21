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
field_axis=[cosd(46.5);sind(46.5);0];
times=(0:60)'*1e-3;
abs_tol=1e-4;
n_configs=50;
pool=gcp('nocreate');
n_workers=0;
if ~isempty(pool)

    % Use the Spinach pool if the caller has started one
    n_workers=pool.NumWorkers;

end

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

% Evaluate the paper's dipolar prefactors in rad s^-1 m^3
mu0=1.25663706127e-6;
hbar=6.62607015e-34/(2*pi);
mu_bohr=9.2740100657e-24;
gamma_w=spin('183W');
g_perp=8.38;
hyperfine=mu0/(4*pi)*g_perp*mu_bohr*gamma_w;
dipolar=mu0/(4*pi)*hbar*gamma_w^2;

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

% Propagate each independently sampled bath before averaging the echoes
realisations=zeros(numel(times),n_configs,2);
parfor (n=1:n_configs,n_workers)

    % Evaluate the two independent bath configurations assigned to this worker
    lattice_echo=cce_secular_pairs(occupied{n},field_axis,hyperfine,dipolar,times,abs_tol);
    random_echo=cce_secular_pairs(amorphous{n},field_axis,hyperfine,dipolar,times,abs_tol);
    realisations(:,n,:)=reshape([lattice_echo random_echo],numel(times),1,2);

end
echo=reshape(mean(realisations,2),numel(times),2);

% Plot the normalised amplitude, not its square
kfigure();
plot(1e3*times,echo,'LineWidth',1.5);
kxlabel('total echo time, $2\tau$ / ms');
kylabel('normalised Hahn echo amplitude');
klegend({'lattice bath','amorphous bath'},'Location','northeast');
kgrid();

end


