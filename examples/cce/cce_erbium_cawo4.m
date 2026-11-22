% Reconstructs the lattice and amorphous Er:CaWO4 Hahn echoes in Fig. S12
% of Le Dantec et al., doi:10.1126/sciadv.abj9786. The curves are
% ensemble means over 50 independently sampled baths of each type.
% Numerical inputs are from the authors' Fortran and supplemental theory
% in doi:10.7910/DVN/9VSAUT. Original bath realisations and seeds were not
% archived. The amorphous sampling in the authors' Fortran sets no minimum
% distance and lets nuclei overlap; this is corrected here by rejecting
% continuum sites closer than 1.12 Angstrom, twice the Shannon crystal
% radius of four-coordinate W(6+).
%
% ilya.kuprov@weizmann.ac.il

function cce_erbium_cawo4()

% Experimental geometry and total Hahn echo times
cell_size=1e-10*[5.2429 5.2429 11.3737];
radius=11e-9;
abundance=0.145;
min_dist=2*0.56e-10;
parameters.observed=1;
parameters.transition=struct('spins',1,'ms',[1/2 -1/2]);
parameters.orientation=deg2rad([0 90 133.5]);
parameters.times=(0:60)*1e-3;
parameters.cce.approximation='projected';
parameters.cce.order=2;
parameters.cce.backend='analytical_hahn';
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

% Effective electron g factor is specified for each physical bath below
sys.tols.prox_cutoff=0;
sys.magnet=0.067;
sys.output='hush';

% Prepare independent isotope occupancies and continuum configurations
occupied=cell(n_configs,1);
amorphous=cell(n_configs,1);
for n=1:n_configs

    % Select isotope sites by Bernoulli occupancy
    stream=RandStream('mrg32k3a','Seed',1);
    stream.Substream=n;
    occupied{n}=lattice(rand(stream,size(lattice,1),1)<abundance,:);

    % Add uniform continuum sites one at a time
    amorphous{n}=zeros(n_amorphous,3);
    k=0;
    while k<n_amorphous

        % Draw a trial site uniformly within the bath sphere
        direction=randn(stream,1,3);
        site=radius*rand(stream)^(1/3)*direction/norm(direction,2);

        % Accept the site if it respects the hard core
        if all(vecnorm(amorphous{n}(1:k,:)-site,2,2)>=min_dist)
            k=k+1;
            amorphous{n}(k,:)=site;
        end

    end

end

% Run every bath realisation through the Hahn echo
baths=[occupied; amorphous];
echoes=zeros(numel(parameters.times),numel(baths));
parfor n=1:numel(baths)

    % Central electron first, then the tungsten bath in Angstrom
    localsys=sys;
    localsys.isotopes=[{'E'} repmat({'183W'},1,size(baths{n},1))];
    inter=struct();
    inter.zeeman.scalar=cell(size(localsys.isotopes));
    inter.zeeman.scalar{1}=8.38;
    inter.coordinates=num2cell([zeros(1,3); 1e10*baths{n}],2);

    % Spinach housekeeping and the echo
    spin_system=create(localsys,inter);
    echoes(:,n)=clexpn(spin_system,@hahn_echo,parameters,'labframe');

end

% Average the lattice and the amorphous ensembles
echo=[mean(echoes(:,1:n_configs),2) mean(echoes(:,(n_configs+1):end),2)];

% Plot the normalised amplitude, not its square
kfigure();
plot(1e3*parameters.times,echo,'LineWidth',1.5);
kxlabel('total echo time, $2\tau$ / ms');
kylabel('normalised Hahn echo amplitude');
klegend({'lattice bath','amorphous bath'},'Location','northeast');
kgrid();

end


