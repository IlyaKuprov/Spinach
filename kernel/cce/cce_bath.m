% Construct explicit, crystalline or random spin baths. Syntax:
%
%                  [sys,inter,bath]=cce_bath(specification)
%
% Required: mode ('explicit','lattice','random'), isotopes (cell row),
% magnetic (tesla). Explicit coordinates and lattice cell vectors are in
% Angstrom; lattice fractional coordinates are rows, cell vectors are rows,
% repeats is a positive integer triple. Random baths use bounds_A (2-by-3),
% count or density (per metre^dimension), dimension, and an explicit seed.
% Occupancy is a site-by-species probability matrix with species names;
% missing probability is a vacancy. Joint configurations instead supply a
% struct array configurations(isotopes,weight). All sampling uses local RNG.
% Optional rotation (proper Cartesian matrix), translation_A, remove_ids,
% selection (centres_A,inner_A,outer_A), groups, tensors_Hz and pairs_Hz
% (records ids,tensor) retain stable site IDs. No dipoles are added here.
% isotopologues.mode='enumerate' or 'sample' produces weighted configurations;
% enumerate requires max_configurations; sample requires count and seed.
%
% ilya.kuprov@weizmann.ac.il

function [sys,inter,bath]=cce_bath(specification)

% Validate the explicit physical boundary
s=specification; grumble(s);
if isfield(s,'configurations')&&strcmp(s.mode,'lattice')&&prod(s.repeats)>1
    error('Joint isotopologues require explicit sites, not implicit lattice replication.');
end
stream=[];
if isfield(s,'seed'), stream=RandStream('mt19937ar','Seed',s.seed); end

% Construct deterministic candidate positions
switch s.mode
    case 'explicit'
        xyz=s.coordinates_A; origin=1:size(xyz,1);
    case 'lattice'
        [i,j,k]=ndgrid(0:s.repeats(1)-1,0:s.repeats(2)-1,0:s.repeats(3)-1);
        shifts=[i(:) j(:) k(:)]; xyz=zeros(0,3); origin=[];
        for n=1:size(shifts,1)
            xyz=[xyz;(s.fractional+shifts(n,:))*s.cell_A]; %#ok<AGROW>
            origin=[origin 1:size(s.fractional,1)]; %#ok<AGROW>
        end
    case 'random'
        span=diff(s.bounds_A,1,1); active=find(span>0);
        measure=prod(span(active)*1e-10);
        if isfield(s,'count'), count=s.count;
        else, count=round(s.density*measure); end
        xyz=zeros(count,3); accepted=0; attempts=0;
        while accepted<count&&attempts<s.max_attempts
            candidate=s.bounds_A(1,:)+rand(stream,1,3).*span;
            attempts=attempts+1;
            if isfield(s,'exclusion_A')&&accepted>0&&...
               any(vecnorm(xyz(1:accepted,:)-candidate,2,2)<s.exclusion_A)
                continue;
            end
            accepted=accepted+1; xyz(accepted,:)=candidate;
        end
        if accepted<count, error('cce_bath:packing','Random packing attempt ceiling reached.'); end
        origin=ones(1,count);
end

% Expand site properties and preserve pre-deletion identities
ids=1:size(xyz,1); species=s.isotopes;
if numel(species)==1, species=repmat(species,1,numel(ids));
elseif numel(species)~=numel(ids), species=species(origin); end
prob=[];
if isfield(s,'occupancy')
    prob=s.occupancy;
    if size(prob,1)==1, prob=repmat(prob,numel(ids),1);
    elseif size(prob,1)~=numel(ids), prob=prob(origin,:); end
end

% Apply explicitly specified frame transformations
rotation=eye(3); translation=zeros(1,3);
if isfield(s,'rotation'), rotation=s.rotation; end
if isfield(s,'translation_A'), translation=s.translation_A; end
xyz=xyz*rotation.'+translation;

% Delete sites and apply multi-centre shell union selection
keep=true(size(ids));
if isfield(s,'remove_ids'), keep=~ismember(ids,s.remove_ids); end
if isfield(s,'selection')
    nearest=inf(numel(ids),1);
    for n=1:size(s.selection.centres_A,1)
        nearest=min(nearest,vecnorm(xyz-s.selection.centres_A(n,:),2,2));
    end
    keep=keep&(nearest.'>=s.selection.inner_A)&(nearest.'<=s.selection.outer_A);
end
ids=ids(keep); xyz=xyz(keep,:); species=species(keep);
if ~isempty(prob), prob=prob(keep,:); end
if size(unique(xyz,'rows'),1)~=size(xyz,1), error('cce_bath:duplicates','Duplicate positions.'); end

% Construct weighted isotopologues before choosing a physical realization
configs=struct('ids',{},'isotopes',{},'weight',{},'multiplicities',{});
if isfield(s,'configurations')
    for n=1:numel(s.configurations)
        names=s.configurations(n).isotopes(keep);
        configs(n)=configuration(ids,names,s.configurations(n).weight);
    end
elseif isfield(s,'isotopologues')
    if isempty(prob), error('Isotopologue enumeration requires occupancy.'); end
    options=s.isotopologues;
    if strcmp(options.mode,'enumerate')
        states={cell(1,0)}; weights=1;
        for n=1:numel(ids)
            p=[prob(n,:) 1-sum(prob(n,:))]; choices=find(p>0);
            if numel(states)*numel(choices)>options.max_configurations
                error('cce_bath:configurations','Enumeration ceiling reached.');
            end
            next={}; next_weights=[];
            for a=1:numel(states)
                for b=choices
                    labels=[s.species {''}];
                    next{end+1}=[states{a} labels(b)]; %#ok<AGROW>
                    next_weights(end+1)=weights(a)*p(b); %#ok<AGROW>
                end
            end
            states=next; weights=next_weights;
        end
    else
        states=cell(1,options.count); weights=ones(1,options.count)/options.count;
        for n=1:options.count, states{n}=draw(prob,s.species,stream); end
    end
    for n=1:numel(states), configs(n)=configuration(ids,states{n},weights(n)); end
end

% Sample occupancy only for a requested single physical realization
if ~isempty(prob)&&isempty(configs), species=draw(prob,s.species,stream); end
if ~isempty(configs)
    chosen=configs(1); active=ismember(ids,chosen.ids); species=chosen.isotopes;
else
    active=~cellfun(@isempty,species); species=species(active);
end
candidate_ids=ids; candidate_xyz=xyz; ids=ids(active); xyz=xyz(active,:);

% Translate the selected realization to native Spinach inputs
sys.isotopes=species; sys.magnet=s.magnetic;
inter.coordinates=mat2cell(xyz,ones(size(xyz,1),1),3).';
if isfield(s,'tensors_Hz')
    inter.zeeman.matrix=cell(1,numel(ids));
    for n=1:numel(ids)
        tensor=s.tensors_Hz(:,:,ids(n));
        inter.zeeman.matrix{n}=rotation*tensor*rotation.';
    end
    inter.zeeman.units='Hz';
end
if isfield(s,'pairs_Hz')
    inter.coupling.matrix=cell(numel(ids));
    for n=1:numel(s.pairs_Hz)
        [present,index]=ismember(s.pairs_Hz(n).ids,ids);
        if all(present)
            tensor=s.pairs_Hz(n).tensor;
            if isscalar(tensor), tensor=tensor*eye(3); end
            inter.coupling.matrix{index(1),index(2)}=rotation*tensor*rotation.';
        end
    end
end

% Return explicit provenance, occupancy and immutable physical group labels
bath.schema_version=1; bath.ids=ids; bath.coordinates_A=xyz;
bath.isotopes=species; bath.configurations=configs;
bath.candidate_ids=candidate_ids; bath.candidate_coordinates_A=candidate_xyz;
bath.rotation=rotation; bath.translation_A=translation;
bath.configuration_selection=1;
if isempty(configs), bath.configuration_selection=[]; end
bath.groups=num2cell(ids); bath.units=struct('coordinates','angstrom','couplings','Hz');
if isfield(s,'groups')
    bath.groups=cellfun(@(g)intersect(g,ids,'stable'),s.groups,'UniformOutput',false);
    bath.groups=bath.groups(~cellfun(@isempty,bath.groups));
    if ~isequal(sort([bath.groups{:}]),ids), error('Groups must partition retained IDs.'); end
end
if strcmp(s.mode,'random')
    bath.measure_SI=measure; bath.dimension=s.dimension;
    bath.realized_density=numel(ids)/measure;
end

end

function names=draw(prob,species,stream)

% Draw each categorical site independently, including explicit vacancies
names=cell(1,size(prob,1)); labels=[species {''}];
for k=1:size(prob,1)
    choice=find(rand(stream)<cumsum([prob(k,:) 1-sum(prob(k,:))]),1);
    names{k}=labels{choice};
end

end

function item=configuration(ids,names,weight)

% Preserve the actual Hilbert dimension of every occupied isotopic species
active=~cellfun(@isempty,names); item.ids=ids(active);
item.isotopes=names(active); item.weight=weight;
item.multiplicities=zeros(size(item.ids));
for k=1:numel(item.ids), [~,item.multiplicities(k)]=spin(item.isotopes{k}); end

end

function grumble(s)

% Require an explicit and finite input schema
if ~isstruct(s)||~isscalar(s)||~all(isfield(s,{'mode','isotopes','magnetic'}))||...
   ~ismember(s.mode,{'explicit','lattice','random'})||~iscellstr(s.isotopes)||...
   isempty(s.isotopes)||~isnumeric(s.magnetic)||~isscalar(s.magnetic)||...
   ~isfinite(s.magnetic)||s.magnetic<0
    error('Invalid bath specification.');
end
for k=1:numel(s.isotopes), spin(s.isotopes{k}); end
if isfield(s,'seed')&&(~isscalar(s.seed)||s.seed<0||s.seed~=fix(s.seed))
    error('Seed must be a nonnegative integer.');
end
switch s.mode
    case 'explicit'
        if ~isfield(s,'coordinates_A')||~isnumeric(s.coordinates_A)||...
           size(s.coordinates_A,2)~=3||~isreal(s.coordinates_A)||any(~isfinite(s.coordinates_A),'all')
            error('Finite Cartesian Angstrom coordinates are required.');
        end
        n=size(s.coordinates_A,1);
    case 'lattice'
        if ~all(isfield(s,{'cell_A','fractional','repeats'}))||...
           ~isequal(size(s.cell_A),[3 3])||any(~isfinite(s.cell_A),'all')||...
           abs(det(s.cell_A))==0||size(s.fractional,2)~=3||...
           any(~isfinite(s.fractional),'all')||any(s.fractional<0|s.fractional>=1,'all')||...
           numel(s.repeats)~=3||any(s.repeats<1|s.repeats~=fix(s.repeats))
            error('Invalid nonorthogonal lattice specification.');
        end
        n=size(s.fractional,1);
    case 'random'
        if ~all(isfield(s,{'bounds_A','dimension','seed','max_attempts'}))||...
           ~isequal(size(s.bounds_A),[2 3])||any(~isfinite(s.bounds_A),'all')||...
           any(diff(s.bounds_A)<0)||sum(diff(s.bounds_A)>0)~=s.dimension||...
           ~ismember(s.dimension,1:3)||s.max_attempts<1||s.max_attempts~=fix(s.max_attempts)||...
           (isfield(s,'count')==isfield(s,'density'))
            error('Random baths need bounds, dimension, seed, attempt ceiling, and count OR density.');
        end
        if isfield(s,'count')&&(~isscalar(s.count)||s.count<0||s.count~=fix(s.count)), error('Invalid count.'); end
        if isfield(s,'density')&&(~isscalar(s.density)||~isfinite(s.density)||s.density<0), error('Invalid SI density.'); end
        if isfield(s,'exclusion_A')&&(~isscalar(s.exclusion_A)||~isfinite(s.exclusion_A)||s.exclusion_A<0), error('Invalid exclusion.'); end
        n=1;
end
if numel(s.isotopes)~=1&&numel(s.isotopes)~=n, error('Species count mismatch.'); end
if isfield(s,'selection')
    q=s.selection;
    if ~isstruct(q)||~all(isfield(q,{'centres_A','inner_A','outer_A'}))||...
       isempty(q.centres_A)||size(q.centres_A,2)~=3||any(~isfinite(q.centres_A),'all')||...
       ~isscalar(q.inner_A)||~isscalar(q.outer_A)||q.inner_A<0||q.outer_A<q.inner_A||...
       isnan(q.inner_A)||isnan(q.outer_A), error('Invalid shell selection.'); end
end
if isfield(s,'remove_ids')&&(~isnumeric(s.remove_ids)||any(~isfinite(s.remove_ids))||...
   any(s.remove_ids<1|s.remove_ids~=fix(s.remove_ids))), error('Invalid deleted site IDs.'); end
if isfield(s,'tensors_Hz')&&(~isnumeric(s.tensors_Hz)||~isreal(s.tensors_Hz)||...
   size(s.tensors_Hz,1)~=3||size(s.tensors_Hz,2)~=3||any(~isfinite(s.tensors_Hz),'all'))
    error('Local tensors must be finite real 3-by-3-by-site arrays in Hz.');
end
if isfield(s,'pairs_Hz')
    for k=1:numel(s.pairs_Hz)
        pair=s.pairs_Hz(k);
        if ~all(isfield(pair,{'ids','tensor'}))||numel(pair.ids)~=2||...
           any(~isfinite(pair.ids)|pair.ids<1|pair.ids~=fix(pair.ids))||...
           ~isnumeric(pair.tensor)||~isreal(pair.tensor)||any(~isfinite(pair.tensor),'all')||...
           (~isscalar(pair.tensor)&&~isequal(size(pair.tensor),[3 3]))
            error('Invalid Hz pair tensor record.');
        end
    end
end
if isfield(s,'rotation')&&(~isequal(size(s.rotation),[3 3])||...
   any(~isfinite(s.rotation),'all')||norm(s.rotation.'*s.rotation-eye(3),'fro')>100*eps||...
   abs(det(s.rotation)-1)>100*eps), error('Rotation must be proper orthogonal.'); end
if isfield(s,'translation_A')&&(~isequal(size(s.translation_A),[1 3])||any(~isfinite(s.translation_A))), error('Invalid translation.'); end
if isfield(s,'occupancy')
    if ~isfield(s,'species')||~iscellstr(s.species)||size(s.occupancy,2)~=numel(s.species)||...
       ~ismember(size(s.occupancy,1),[1 n])||any(~isfinite(s.occupancy),'all')||...
       any(s.occupancy<0,'all')||any(sum(s.occupancy,2)>1)
        error('Invalid isotope probabilities.');
    end
    for k=1:numel(s.species), spin(s.species{k}); end
    if ~isfield(s,'isotopologues')&&~isfield(s,'seed'), error('Occupancy sampling requires a seed.'); end
end
if isfield(s,'configurations')
    weights=[s.configurations.weight];
    if isempty(weights)||any(~isfinite(weights)|weights<0)||abs(sum(weights)-1)>numel(weights)*eps
        error('Joint configuration weights must sum to one.');
    end
    for k=1:numel(weights)
        if numel(s.configurations(k).isotopes)~=n, error('Joint configuration size mismatch.'); end
    end
end
if isfield(s,'isotopologues')
    o=s.isotopologues;
    if ~isfield(o,'mode')||~ismember(o.mode,{'enumerate','sample'}), error('Invalid isotopologue mode.'); end
    if strcmp(o.mode,'enumerate')
        if ~isfield(o,'max_configurations')||o.max_configurations<1||o.max_configurations~=fix(o.max_configurations), error('Enumeration ceiling required.'); end
    elseif ~isfield(s,'seed')||~isfield(o,'count')||o.count<1||o.count~=fix(o.count)
        error('Sampling count and seed required.');
    end
end

end
