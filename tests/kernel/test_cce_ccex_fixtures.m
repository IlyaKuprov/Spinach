% Pinned CCEX analytic and reduced-basic01 reference comparisons. Syntax:
%
%                     result=test_cce_ccex_fixtures()
%
% CCEX cbcde24292e0e99755967267fe7c477fc029cb4d, CPU executable.
% Saved decimal traces print ten digits after the point. The basic01
% reduction uses the nearest six sites from its 978-spin selected bath.
% Input units are Angstrom, Gauss, ms and rad/(ms Gauss); no isotope
% constants are silently substituted for the reference's declared values.
%
% ilya.kuprov@weizmann.ac.il

function result=test_cce_ccex_fixtures()

% Load externally executed complex reference trajectories
result=new_test_result('kernel/cce_ccex_fixtures','CCEX fixture parity',...
                       'Analytic FID/Hahn and basic01 nearest-six CCE2.');
data=jsondecode(fileread(fullfile(fileparts(mfilename('fullpath')),...
                                'test_cce_ccex_data.json')));
errors=zeros(1,numel(data));
for f=1:numel(data)
    if iscell(data), fixture=data{f}; else, fixture=data(f); end
    xyz=fixture.xyz; if isvector(xyz), xyz=reshape(xyz,1,3); end
    nb=size(xyz,1); center=fixture.center(:).';

    % Construct the source's full Cartesian point-dipole tensors
    sys.isotopes=[{'14N'} repmat({'13C'},1,nb)]; sys.magnet=0;
    sys.output='hush'; sys.disable={'hygiene'}; sys.parallel={'local',1};
    inter.zeeman.matrix=repmat({zeros(3)},1,nb+1); inter.coupling.matrix=cell(nb+1);
    inter.coordinates=mat2cell([center;xyz],ones(1,nb+1),3).';
    spin_system=create(sys,inter);
    for j=1:nb+1
        for k=1:nb+1, spin_system.inter.coupling.matrix{j,k}=zeros(3); end
    end
    gamma_n=6.728284; gamma_e=-17608.597050;
    for j=1:nb
        spin_system.inter.zeeman.matrix{j+1}(3,3)=-gamma_n*500*1000;
        r=xyz(j,:)-center;
        spin_system.inter.coupling.matrix{1,j+1}=dipolar(r,gamma_e,gamma_n);
        for k=j+1:nb
            r=xyz(j,:)-xyz(k,:);
            spin_system.inter.coupling.matrix{j+1,k+1}=dipolar(r,gamma_n,gamma_n);
        end
    end

    % Use ms=+1 and ms=0 branches and a selective logical pi pulse
    p.central_spins=1; p.bath_spins=2:nb+1; p.method='cce'; p.expansion='standard';
    p.quantity='coherence'; p.order=fixture.order; p.central.alpha=[1;0;0];
    p.central.beta=[0;1;0]; p.assumptions='labframe'; p.orientation=[0 0 0];
    p.interaction_cutoff=Inf; p.graph=struct('mode','distance','radius',fixture.rdip,...
                                         'score','coupling_norm');
    p.meanfield='none'; p.times=fixture.times_ms(:).'/1000;
    p.bath_state=struct('mode','mixed_product'); p.sequence=[];
    p.numerics=struct('atol',1e-7,'rtol',1e-7,'max_memory_bytes',2^29);
    if fixture.pulses==1
        a=p.central.alpha; b=p.central.beta;
        p.sequence=struct('type','ideal','target',1,'placement','fraction','at',0.5,...
                          'operator',pi*(a*b'+b*a')/2);
    end
    actual=cce(spin_system,p); reference=fixture.real+1i*fixture.imag;
    errors(f)=max(abs(actual.value-reference));
    if isfield(fixture,'raw')
        for j=1:numel(fixture.raw)
            expected=fixture.raw(j); sites=sort(expected.sites(:).');
            index=find(arrayfun(@(r) isequal(r.sites,sites),actual.raw{1}),1);
            record=actual.raw{1}(index); ref=expected.real+1i*expected.imag;
            result=test_true(result,'CCEX raw subset',...
                max(abs(record.value-ref))<=max(record.abs_error_bound)+sqrt(2)*0.5e-10,...
                'Raw singleton and pair traces agree before expansion.');
        end
    end
    budget=max(actual.abs_error_bound)+sqrt(2)*0.5e-10;
    result=test_true(result,'CCEX printed curve',errors(f)<=budget,...
                     'Half the last printed unit plus measured propagation/reduction bound.');
end
result.metrics.maximum_errors=errors;

end

% Reproduce the reference's published point-dipole coefficient in rad/s
function tensor=dipolar(r,g1,g2)
distance=norm(r);
tensor=1000*1.054571729*g1*g2*(eye(3)-3*(r.'*r)/distance^2)/distance^3;
end

