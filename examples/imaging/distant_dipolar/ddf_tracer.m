% Heteronuclear tracer in the dipolar field of a proton grating. The
% synthetic initial proton longitudinal density varies along Z; the
% carbon polarisation is uniform. Both species have transverse moments
% in their own frames. Compare the raw phases at the high-density end
% of the grating for passive and magnetically active carbon. Transport
% reflects at box faces while the field is periodic with needle mean.

function ddf_tracer()

% Specify two uncoupled molecular species with equal reference populations
sys.magnet=1;
sys.isotopes={'1H','13C'};
sys.parallel={'processes',1};
inter.chem.parts={1,2};
inter.chem.concs=[1 1];
bas.formalism='sphten-liouv';
bas.approximation={'none','none'};
spin_system=basis(create(sys,inter),bas);

% Set the two rotating frames and the Cartesian grid
p.spins=sys.isotopes;
p.offset=[0 0];
p.groups={1,2};
p.frames=[1 2];
p.npts=[10 10 10];
p.dims=[1 1 0.15]*1e-3;
p.deriv={'period',3};
p.diff=0;
spacing=p.dims./p.npts;
symbol=ddf_periodic(p.npts,spacing,0);
p.field=@(m)ddf_cartesian(m,p.npts,symbol);

% Prepare the grating and detect one face with a uniform receiver phase
grating=1+0.5*cos(pi*((1:p.npts(3))-0.5)/p.npts(3));
p.rho0_ph={repmat(reshape(grating,1,1,[]),[p.npts(1:2) 1]),ones(p.npts)};
receiver=zeros(p.npts);
receiver(:,:,1)=1/prod(p.npts(1:2));
p.coil_ph={receiver};
dim=spin_system.bas.offsets(end);
p.rlx_ph={ones(p.npts)};
p.rlx_op={sparse(dim,dim)};

% Assign separate reflecting diffusion rates to the two substance blocks
diffusivity=[2e-9 7e-10];
p.transport=sparse(dim*prod(p.npts),dim*prod(p.npts));
ref=zeros(1,2);
for q=1:2
    sites=(spin_system.bas.offsets(q)+1):spin_system.bas.offsets(q+1);
    proj=sparse(sites,sites,1,dim,dim);
    rate=ddf_diffusion(true(p.npts),spacing,diffusivity(q));
    p.transport=p.transport+1i*kron(rate,proj);
    z=coil_state(spin_system,'Lz',q);
    ref(q)=real(z'*z);
    p.rho0_st{q}=z+coil_state(spin_system,'Lx',q);
end

% Sample the free induction at fixed receiver phase every 5 ms
p.duration=[0;0.005*ones(8,1)];
p.acquire=true(9,1);
p.steps=2*ones(9,1);
p.rf=zeros(9,2);
p.gradient=zeros(9,3);

% Compare reciprocal phase accumulation on each receiver independently
kfigure();
for q=1:2
    p.coil_st={coil_state(spin_system,'L+',q)};
    subplot(1,2,q);
    hold on;
    for carbon=[0 0.03]
        p.scale=[0.05 carbon]./ref;
        fid=imaging(spin_system,@ddf_sequence,p);
        plot(1e3*cumsum(p.duration),unwrap(angle(fid)));
    end
    kxlabel('time, ms');
    kylabel([sys.isotopes{q} ' phase, rad']);
    klegend({'passive carbon','active carbon'});
    kgrid();
end

end


