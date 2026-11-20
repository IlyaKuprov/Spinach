% Unrestricted diffusion control for finite rectangular gradient lobes.
% The ratio to a separately simulated zero-gradient acquisition must
% reproduce Stejskal and Tanner, https://doi.org/10.1063/1.1695690 .
% Syntax:
%
%                       result=test_ddf_pgse()
%
% Output:
%
%   result - complex echo ratio comparisons for two encoding axes

function result=test_ddf_pgse()

% Set up a single proton with Cartesian axes defined in the laboratory
result=new_test_result('kernel/ddf_pgse','Finite gradient diffusion',...
                       'Test evolution during rectangular lobes with feedback off.');
sys.magnet=9.4;
sys.isotopes={'1H'};
sys.parallel={'processes',1};
sys.output='hush';
sys.disable={'hygiene'};
bas.formalism='sphten-liouv';
bas.approximation={'none'};
s=basis(create(sys,struct()),bas);

% Barros doped-water diffusivity with diagnostic rectangular pulse timings
diffusivity=2.15e-9;
strength=0.12;
delta=0.004;
separation=0.010;
wave_number=spin('1H')*strength*delta;
target=exp(-diffusivity*wave_number^2*(separation-delta/3));

% Prepare a transverse state and a balanced pair without RF refocusing
p.groups={1};
p.frames=1;
p.scale=0;
p.field=@(magn)zeros(size(magn));
p.duration=[delta;separation-delta;delta];
p.steps=ones(3,1);
p.acquire=[false;false;true];
p.rf=zeros(3,1);
p.spins={'1H'};
p.offset=0;
p.diff=0;
p.deriv={'period',3};
p.rho0_st={state(s,'Lx','1H')};
p.coil_st={coil_state(s,'L+','1H')};
p.rlx_ph={};
p.rlx_op={};

% Rotate the slab and receiver together from transverse to axial encoding
for axis_id=[1 3]
    p.npts=[10 10 10];
    p.npts(axis_id)=60;
    spacing=2e-6;
    p.dims=spacing*p.npts;
    mask=true(p.npts);
    p.rho0_ph={double(mask)};
    rate=ddf_diffusion(mask,p.dims./p.npts,diffusivity);
    p.transport=1i*kron(rate,speye(s.bas.offsets(end)));

    % Detect the interior, five diffusion lengths from the reflecting walls
    centre={':',':',':'};
    centre{axis_id}=21:40;
    receiver=zeros(p.npts);
    receiver(centre{:})=1;
    p.coil_ph={receiver/sum(receiver(:))};

    % Keep the full acquisition time and transport in the gain reference
    p.gradient=zeros(3,3);
    reference=imaging(s,@ddf_sequence,p);
    p.gradient(:,axis_id)=[strength;0;-strength];
    echo=imaging(s,@ddf_sequence,p);
    ratio=echo/reference;

    % Bound spatial discretisation error below the instantaneous-kick bias
    label=sprintf('axis %d, echo %.6f, theory %.6f',axis_id,real(ratio),target);
    result=test_close(result,label,ratio,target,0.002,0,...
                      'finite rectangular gradient moments and in-lobe diffusion');

    % Integrate the face-flux lattice symbol along the gradient trajectory
    phase=wave_number*spacing;
    area=2*delta*(sin(phase)/phase-1)+(separation-delta)*(cos(phase)-1);
    lattice=exp(2*diffusivity*area/spacing^2);
    result=test_close(result,sprintf('axis %d, lattice oracle',axis_id),...
                      ratio,lattice,2e-6,0,'separate spatial error from propagation error');
end

end


