% Projected CCE-1 Hahn echo of the axial (kk) divacancy in natural 4H-SiC.
% Reconstructs Supplementary Figure 5 of Seo et al. (10.1038/ncomms12935)
% with a point electron dipole at the Si vacancy and independent isotope
% occupancies rather than the authors' finite ensemble of random baths.
% Nuclear pair correlations and the associated echo decay are absent.
% Calculation time: about 12 minutes on a Xeon Platinum 8270 (serial CCE).
%
%             [answer,diagnostics]=cce_sic_divacancy()
%
% answer is the normalised echo (time rows, magnetic field columns).
% diagnostics contains the time grid (s), fields (T), site count, and
% number of distinct occupied-site calculations per field.
%
% ilya.kuprov@weizmann.ac.il

function [answer,diagnostics]=cce_sic_divacancy()

% Hexagonal 4H-SiC cell in Angstrom and fractional Si/C basis
lattice_a=3.073;
lattice_c=10.053;
sites=[0   0   0.188;
       0   0   0.688;
       1/3 2/3 0.438;
       2/3 1/3 0.938;
       0   0   0.000;
       0   0   0.500;
       1/3 2/3 0.250;
       2/3 1/3 0.750];
isotopes={'29Si','13C'};
abundance=[0.047 0.011];

% Generate cell translations covering the 50 Angstrom bath
radius=50;
nxy=ceil(2*radius/(sqrt(3)*lattice_a));
nz=ceil(radius/lattice_c);
[tx,ty,tz]=ndgrid(-nxy:nxy,-nxy:nxy,-nz:nz);
translations=[tx(:) ty(:) tz(:)];
fractional=repelem(translations,8,1)+repmat(sites,size(translations,1),1);
species=repmat([ones(4,1);2*ones(4,1)],size(translations,1),1);

% Remove both axial vacancy sites and centre the dipole on the silicon vacancy
keep=~ismember(fractional,sites([1 5],:),'rows');
fractional=fractional(keep,:)-sites(1,:);
species=species(keep);

% Group equal singleton responses by isotope, basal distance, and absolute height
basal=round(3*fractional(:,1:2));
height=abs(round(1000*fractional(:,3)));
basal_sq=basal(:,1).^2+basal(:,2).^2-basal(:,1).*basal(:,2);
rho=lattice_a*sqrt(basal_sq)/3;
z=lattice_c*height/1000;
keep=hypot(rho,z)<=radius;
[~,representatives,group]=unique([species(keep) basal_sq(keep) height(keep)],'rows');
xyz=[rho(keep) zeros(nnz(keep),1) z(keep)];
xyz=xyz(representatives,:);
species=species(keep);
species=species(representatives);
counts=accumarray(group,1);

% Use spin-one electron manifolds and maximally mixed nuclear states
sys.output='hush';
sys.parallel={'processes',1};
bas.formalism='zeeman-hilb';
bas.approximation='none';
parameters.observed=1;
parameters.transition.spins=1;
parameters.transition.ms=[1 0];
parameters.times=0:1e-6:2e-3;
parameters.cce.order=1;
parameters.cce.approximation='projected';
parameters.cce.graph=struct('mode','complete');
parameters.cce.bath_state.mode='mixed_product';
parameters.cce.basis=bas;
parameters.cce.numerics=struct('atol',1e-8,'rtol',1e-8,...
                               'max_memory_bytes',2^28);

% Evaluate independent sites in batches to bound the dense pair tables
fields=[25 65 340]*1e-4;
batch_size=128;
answer=ones(numel(parameters.times),numel(fields));
for k=1:numel(fields)

    % Set the c-axis field for this trace
    sys.magnet=fields(k);
    for n=1:batch_size:numel(counts)

        % Build the projected electron-secular couplings in Hz
        batch=n:min(n+batch_size-1,numel(counts));
        sys.isotopes=[{'E3'} isotopes(species(batch))];
        inter.coupling.matrix=cell(numel(batch)+1);
        for j=1:numel(batch)

            % Convert the point-dipole tensor and retain only its electron-z row
            site=batch(j);
            tensor=xyz2hfc([0 0 0],xyz(site,:),isotopes{species(site)});
            inter.coupling.matrix{1,j+1}=[zeros(2,3);...
                                         -spin('E')*1e-4*tensor(3,:)/(2*pi)].';

        end

        % Get the native projected Hahn echo for every occupied site
        parameters.cce.bath_spins=2:numel(sys.isotopes);
        spin_system=create(sys,inter);
        [~,details]=clexpn(spin_system,@hahn_echo,parameters,'labframe');

        % Average each independent occupied-or-empty factor before multiplying
        occupied=[details.raw{1}(2:end).value];
        weights=abundance(species(batch));
        factors=1+(occupied-1).*weights;
        multiplicity=reshape(counts(batch),1,[]);
        answer(:,k)=answer(:,k).*prod(factors.^multiplicity,2);

    end

end

% Return experimental axes and crystallographic site counts
diagnostics=struct('times',parameters.times(:),'fields',fields,...
                   'site_count',sum(counts),'distinct_sites',numel(counts));

% Offset the traces as in the published supplementary figure
figure();
plot(1e3*diagnostics.times,real(answer)+[2 1 0]);
xlabel('total echo time / ms');
ylabel('coherence (traces offset)');
legend('25 G','65 G','340 G','Location','eastoutside');
kgrid();

end


