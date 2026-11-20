% Independent nonlinear Bloch reference for overlapping controls. Syntax:
%
%                        result=test_ddf_waveform()
%
% Outputs:
%
%   result - complex samples at three finite waveform edges

function result=test_ddf_waveform()

% Build a single-spin system and use a phase-sensitive local receiver
result=new_test_result('kernel/ddf_waveform','Finite DDF waveform',...
                       'Compare concurrent RF, gradient, and dipolar evolution with Bloch integration.');
sys.magnet=1;
sys.isotopes={'1H'};
sys.parallel={'processes',1};
sys.output='hush';
sys.disable={'hygiene'};
bas.formalism='sphten-liouv';
bas.approximation={'none'};
s=basis(create(sys,struct()),bas);
p.npts=[10 10 10];
p.dims=[1 2 3]*1e-3;
p.deriv={'period',3};
p.diff=0;
p.spins={'1H'};
p.offset=0;
p.groups={1};
p.frames=1;
z=coil_state(s,'Lz',1);
p.rho0_st={z};
p.rho0_ph={ones(p.npts)};
p.coil_st={coil_state(s,'L+',1)};
p.coil_ph={zeros(p.npts)};
p.coil_ph{1}(1,1,8)=1;
p.rlx_ph={ones(p.npts)};
p.rlx_op={sparse(numel(z),numel(z))};
dim=prod(p.npts)*numel(z);
p.transport=sparse(dim,dim);

% Match the mean to the axial-wave symbol so that feedback is local
symbol=ddf_periodic(p.npts,p.dims./p.npts,1);
p.field=@(m)ddf_cartesian(m,p.npts,symbol);
density=0.04;
p.scale=density/real(z'*z);
p.duration=[2;3;5]*1e-3;
p.steps=[12;18;30];
p.acquire=true(3,1);
p.gradient=[0 0 1;0 0 -0.5;0 0 0]*1e-3;
p.rf=[300+200i;-100+40i;0];
fid=imaging(s,@ddf_sequence,p);

% Integrate the torque of the local dipolar field at the selected voxel
gamma=spin('1H');
coord=(8-(p.npts(3)+1)/2)*p.dims(3)/p.npts(3);
dipolar=gamma*1.25663706127e-6*density;
moment=[0;0;1];
expected=zeros(3,1);
for n=1:3
    external=[real(p.rf(n));imag(p.rf(n));-gamma*coord*p.gradient(n,3)];
    rhs=@(~,m)cross(external+[0;0;dipolar*m(3)],m);
    [~,trajectory]=ode45(rhs,[0 p.duration(n)],moment,...
                         odeset('AbsTol',1e-12,'RelTol',1e-10));
    moment=trajectory(end,:).';
    expected(n)=real(z'*z)*(moment(1)-1i*moment(2));
end
result=test_close(result,'three waveform edges',fid,expected,2e-5,0,...
                  'the axial tensor reduces to an Mz-dependent offset for any tilt');

% Isolate resonant uniform-coil radiation damping from the dipolar field
p.rho0_st={coil_state(s,'Lx',1)};
p.gradient=zeros(1,3);
p.rf=0;
p.acquire=true;
p.duration=0.04;
p.steps=160;
damping=1/(gamma*density*p.duration);
p.field=@(m)repmat(damping*[-mean(m(:,2)) mean(m(:,1)) 0],size(m,1),1);
damped=imaging(s,@ddf_sequence,p);
p.scale=0;
control=imaging(s,@ddf_sequence,p);
result=test_close(result,'uniform-coil radiation damping',damped/control,sech(1),...
                  1e-6,0,'an initially transverse moment relaxes as sech(t/T_rd)');

end


