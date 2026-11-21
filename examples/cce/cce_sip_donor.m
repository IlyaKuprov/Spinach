% Orientation-dependent Hahn echo of a phosphorus donor in natural silicon,
% following Onizhuk and Galli, doi:10.1002/adts.202100254, Fig. 5(a).
% A fixed isotope realisation is used with secular CCE2 and an unpolarised
% bath. The field and pair cutoff follow the archived PyCCE tutorial; the
% original Fig. 5(a) isotope configuration was not published.
%
% Syntax:
%
%            [answer,diagnostics]=cce_sip_donor()
%
% Returns a complex coherence matrix (time by orientation) on a uniform
% 0--1 ms grid, with cluster counts and validity for each orientation.
%
% Calculation time: 7 minutes (six Xeon Platinum 8270 workers).
%
% ilya.kuprov@weizmann.ac.il

function [answer,diagnostics]=cce_sip_donor()

% Natural silicon diamond lattice with P on the origin site
lattice=5.43;
bath_radius=80;
pair_radius=10;
extent=ceil(bath_radius/lattice);
[gx,gy,gz]=ndgrid(-extent:extent);
fcc=[0 0 0;0 2 2;2 0 2;2 2 0]/4;
sites=[fcc;fcc+1/4];
xyz=lattice*(repmat([gx(:) gy(:) gz(:)],8,1)+...
             repelem(sites,numel(gx),1));
radius=vecnorm(xyz,2,2);
xyz=xyz((radius>0)&(radius<=bath_radius),:);
stream=RandStream('mt19937ar','Seed',8800);
xyz=xyz(rand(stream,size(xyz,1),1)<0.047,:);
radius=vecnorm(xyz,2,2);
count=size(xyz,1);

% Donor contact coupling in Hz from the six Kohn--Luttinger valleys
long_axis=0.81*25.09;
short_axis=0.81*14.43;
wavevector=0.85*2*pi/lattice;
envelope=exp(-sqrt(xyz.^2/short_axis^2+...
                   (radius.^2-xyz.^2)/long_axis^2));
amplitude=sum(envelope.*cos(wavevector*xyz),2)/...
          sqrt(pi*long_axis^2*short_axis);
dip_scale=1e-7*1.054571817e-34*1e30/(2*pi);
gamma_e=spin('E');
gamma_si=spin('29Si');
contact=(16*pi/9)*186*dip_scale*gamma_e*gamma_si*amplitude.^2;
dip_mask=radius>=long_axis;

% Retain bath pairs within the distance cutoff
distance=pdist2(xyz,xyz);
[pair_a,pair_b]=find(triu(distance<=pair_radius,1));
displacement=xyz(pair_a,:)-xyz(pair_b,:);
pair_dist=vecnorm(displacement,2,2);
pair_unit=displacement./pair_dist;
pair_scale=dip_scale*gamma_si^2./pair_dist.^3;

% Use the rotating frame of both Zeeman terms after the high-field projection
sys.magnet=0;
sys.isotopes={'E'};
sys.output='hush';
sys.parallel={'processes',6};
parameters.observed=1;
parameters.transition.spins=1;
parameters.transition.ms=[1/2 -1/2];
parameters.times=linspace(0,1e-3,101);
parameters.cce.order=2;
parameters.cce.approximation='projected';
parameters.cce.bath_state.mode='mixed_product';
parameters.cce.basis=struct('formalism','zeeman-hilb','approximation','none');
parameters.cce.numerics=struct('atol',1e-5,'rtol',0,...
                               'max_memory_bytes',2^28);

% Tile the distance-selected pairs without splitting any pair
pairs=[pair_a pair_b];
isolated=setdiff((1:count).',pairs(:));
tiles=[pairs;isolated isolated];
tile_size=128;
tile_count=ceil(size(tiles,1)/tile_size);
tile_ids=cell(1,tile_count);
tile_rows=cell(1,tile_count);
duplicate=cell(1,tile_count);
seen=false(count,1);
for k=1:tile_count

    % Assign each physical monomer to its first tile
    rows=(k-1)*tile_size+1:min(k*tile_size,size(tiles,1));
    ids=unique(tiles(rows,:));
    tile_ids{k}=ids;
    tile_rows{k}=rows(rows<=size(pairs,1));
    duplicate{k}=find(seen(ids))+1;
    seen(ids)=true;

end

% Initialise the pool once before dispatching independent cluster tiles
create(sys,struct());
directions=[0 0 1;1 1 1;1 1 0];
answer=zeros(numel(parameters.times),size(directions,1));
diagnostics=cell(1,size(directions,1));
for n=1:size(directions,1)

    % Rotate from crystal axes to the field frame
    field=directions(n,:)/norm(directions(n,:));
    field_cos=xyz*field.'./radius;
    pair_cos=pair_unit*field.';

    % Retain the contact and dipolar longitudinal hyperfine terms
    hyperfine=contact+dip_scale*gamma_e*gamma_si*dip_mask.*...
                      (1-3*field_cos.^2)./radius.^3;
    pair_zz=pair_scale.*(1-3*pair_cos.^2);

    % Pack each exact cluster tile for worker-local evaluation
    tile_input=cell(1,tile_count);
    for k=1:tile_count

        % Map the retained edges from physical bath labels to tile labels
        ids=tile_ids{k};
        rows=tile_rows{k};
        [~,local_ids]=ismember(pairs(rows,:),ids);
        sites_in_tile=numel(ids);
        adjacency=sparse(local_ids(:,1),local_ids(:,2),true,sites_in_tile,sites_in_tile);
        tile_params=parameters;
        tile_params.cce.graph=struct('mode','explicit','adjacency',adjacency|adjacency.');
        tile_sys=sys;
        tile_sys.isotopes=[{'E'} repmat({'29Si'},1,sites_in_tile)];
        inter.coupling.matrix=cell(sites_in_tile+1);

        % Pure-dephasing electron coupling to longitudinal nuclear moments
        for j=1:sites_in_tile

            % Secular electron–nuclear tensor in the magnetic field frame
            inter.coupling.matrix{1,j+1}=diag([0 0 hyperfine(ids(j))]);

        end

        % Install only this tile's nuclear dipolar edges
        for j=1:numel(rows)

            % Secular homonuclear dipole includes the nuclear flip-flop term
            dzz=pair_zz(rows(j));
            tensor=diag([-dzz/2 -dzz/2 dzz]);
            inter.coupling.matrix{local_ids(j,1)+1,local_ids(j,2)+1}=tensor;

        end
        tile_input{k}=struct('sys',tile_sys,'inter',inter,'params',tile_params,...
                             'duplicate',duplicate{k});

    end

    % Distribute the native CCE2 tiles with exact monomer de-duplication
    factors=zeros(numel(parameters.times),tile_count);
    valid=false(numel(parameters.times),tile_count);
    parfor k=1:tile_count

        % Evaluate the native CCE2 and divide out already counted monomers
        input=tile_input{k};
        spin_system=create(input.sys,input.inter);
        [tile_echo,tile_diag]=clexpn(spin_system,@hahn_echo,input.params,'labframe');
        repeated=arrayfun(@(c) c.order==1&&ismember(c.sites,input.duplicate),...
                          tile_diag.clusters);
        raw=tile_diag.raw{1};
        raw_monomers=raw(ismember([raw.cluster_id],[tile_diag.clusters(repeated).id]));
        normaliser=prod([ones(numel(tile_echo),1) raw_monomers.value],2);
        factors(:,k)=tile_echo./normaliser;
        valid(:,k)=tile_diag.reduction.valid&isfinite(factors(:,k));

    end

    % The bare ideal Hahn echo is one, so tile factors multiply directly
    answer(:,n)=prod(factors,2);
    diagnostics{n}=struct('times',parameters.times,'valid',all(valid,2),...
                          'bath_count',count,'pair_count',size(pairs,1),...
                          'tile_count',tile_count);

end

% Compare the three normalised Hahn echo amplitudes
kfigure();
plot(1e3*parameters.times,real(answer),'LineWidth',1.5);
kxlabel('Total echo time (ms)');
kylabel('Normalised coherence');
klegend({'[001]','[111]','[110]'},'Location','northeast');
kgrid();

end


