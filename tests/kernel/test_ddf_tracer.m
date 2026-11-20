% Heteronuclear tracer, reciprocal feedback, and transport tests. Syntax:
%
%                         result=test_ddf_tracer()
%
% Outputs:
%
%   result - signed proton/carbon precession and diffusion comparisons

function result=test_ddf_tracer()

% Build independently moving proton and carbon spin-half substances
result=new_test_result('kernel/ddf_tracer','Two-isotope DDF',...
                       'Check source concentration, signed target phase, and species diffusion.');
sys.magnet=1;
sys.isotopes={'1H','13C'};
sys.parallel={'processes',1};
sys.output='hush';
sys.disable={'hygiene'};
inter.chem.parts={1,2};
inter.chem.concs=[1 1];
bas.formalism='sphten-liouv';
bas.approximation={'none','none'};
s=basis(create(sys,inter),bas);

% Prepare both species at 45 degrees in a uniform periodic needle
p.npts=[10 10 10];
p.dims=[1 2 3]*1e-3;
p.deriv={'period',3};
p.diff=0;
p.spins=sys.isotopes;
p.offset=[0 0];
p.groups={1,2};
p.frames=[1 2];
nvox=prod(p.npts);
dim=s.bas.offsets(end);
p.transport=sparse(dim*nvox,dim*nvox);
p.rlx_ph={ones(p.npts)};
p.rlx_op={sparse(dim,dim)};
p.rho0_st=cell(1,2);
p.rho0_ph={ones(p.npts),ones(p.npts)};
ref=zeros(1,2);
for q=1:2
    z=coil_state(s,'Lz',q);
    ref(q)=real(z'*z);
    p.rho0_st{q}=z+coil_state(s,'Lx',q);
end
symbol=ddf_periodic(p.npts,p.dims./p.npts,0);
p.field=@(m)ddf_cartesian(m,p.npts,symbol);
p.coil_ph={ones(p.npts)/nvox};
p.duration=[0;0.01];
p.steps=[1;16];
p.gradient=zeros(2,3);
p.rf=zeros(2,2);
p.acquire=true(2,1);

% Test passive and reciprocal limits with independent fixed receivers
mu0=1.25663706127e-6;
for carbon=[0 0.03]
    density=[0.05 carbon];
    p.scale=density./ref;
    for q=1:2
        p.coil_st={coil_state(s,'L+',q)};
        fid=imaging(s,@ddf_sequence,p);
        rate=spin(sys.isotopes{q})*mu0*(sum(density)+density(q)/2)/3;
        result=test_close(result,sprintf('group %d, carbon %g A/m',q,carbon),...
                          fid(2)/fid(1),exp(1i*rate*p.duration(2)),2e-7,0,...
                          'only the matching frame supplies transverse feedback');
    end
end

% Decay the first reflecting-slab eigenmode at two independent rates
p.scale=[0 0];
diffusivity=[2e-9 7e-10];
spacing=p.dims./p.npts;
profile=repmat(cos(pi*((1:p.npts(1))'-0.5)/p.npts(1)),[1 p.npts(2:3)]);
p.rho0_ph={profile,profile};
p.coil_ph={profile};
p.duration(2)=1;
for q=1:2
    sites=(s.bas.offsets(q)+1):s.bas.offsets(q+1);
    proj=sparse(sites,sites,1,dim,dim);
    rate=ddf_diffusion(true(p.npts),spacing,diffusivity(q));
    p.transport=p.transport+1i*kron(rate,proj);
end
for q=1:2
    p.coil_st={coil_state(s,'L+',q)};
    fid=imaging(s,@ddf_sequence,p);
    rate=diffusivity(q)*(2*sin(pi/(2*p.npts(1)))/spacing(1))^2;
    result=test_close(result,sprintf('diffusion of group %d',q),...
                      fid(2)/fid(1),exp(-rate*p.duration(2)),2e-8,0,...
                      'each molecular block follows its own no-flux Laplacian');
end

% Detect carbon at rest in the field of a diffusing longitudinal proton mode
sites=1:s.bas.offsets(2);
proj=sparse(sites,sites,1,dim,dim);
rate=ddf_diffusion(true(p.npts),spacing,diffusivity(1));
p.transport=1i*kron(rate,proj);
p.scale=[0.05/ref(1) 0];
p.rho0_st={coil_state(s,'Lz',1),coil_state(s,'Lx',2)};
p.rho0_ph={profile,ones(p.npts)};
p.coil_st={coil_state(s,'L+',2)};
p.coil_ph={zeros(p.npts)};
p.coil_ph{1}(1,1,1)=1;
p.duration(2)=0.2;
fid=imaging(s,@ddf_sequence,p);
decay=diffusivity(1)*(2*sin(pi/(2*p.npts(1)))/spacing(1))^2;
phase=spin('13C')*mu0*(0.05/3)*profile(1)*(-expm1(-decay*p.duration(2)))/decay;
result=test_close(result,'tracer in a decaying proton grating',fid(2)/fid(1),...
                  exp(1i*phase),2e-7,0,'reporter phase integrates the evolving source field');

% Rotate only carbon, then encode its phase along all three axes
p.scale=[0 0];
p.transport=sparse(dim*nvox,dim*nvox);
p.rho0_st={coil_state(s,'Lz',1),coil_state(s,'Lz',2)};
p.rho0_ph={ones(p.npts),ones(p.npts)};
p.duration=[1e-3;2e-3];
p.rf=[0 1i*pi/(2*p.duration(1));0 0];
p.gradient(2,:)=[2 -3 5]*1e-3;
p.coil_st={coil_state(s,'L+',2)};
p.coil_ph={zeros(p.npts)};
voxel=[2 4 9];
p.coil_ph{1}(voxel(1),voxel(2),voxel(3))=1;
fid=imaging(s,@ddf_sequence,p);
coord=(voxel-(1+p.npts)/2).*spacing;
phase=spin('13C')*dot(coord,p.gradient(2,:))*p.duration(2);
result=test_close(result,'finite RF and oblique encoding',fid,...
                  ref(2)*[1;exp(1i*phase)],2e-8,0,...
                  'selective Y rotation then cell-centred carbon Zeeman phase');

end


