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
parameters.npts=[10 10 10];
parameters.dims=[1 2 3]*1e-3;
parameters.deriv={'period',3};
parameters.diff=0;
parameters.spins={'1H'};
parameters.offset=0;
parameters.groups={1};
parameters.frames=1;
z=coil_state(s,'Lz',1);
parameters.rho0_st={z};
parameters.rho0_ph={ones(parameters.npts)};
parameters.coil_st={coil_state(s,'L+',1)};
parameters.coil_ph={zeros(parameters.npts)};
parameters.coil_ph{1}(1,1,8)=1;
parameters.rlx_ph={ones(parameters.npts)};
parameters.rlx_op={sparse(numel(z),numel(z))};
dim=prod(parameters.npts)*numel(z);
parameters.transport=sparse(dim,dim);

% Match the mean to the axial-wave symbol so that feedback is local
symbol=ddf_periodic(parameters.npts,parameters.dims./parameters.npts,1);
parameters.field=@(m)ddf_cartesian(m,parameters.npts,symbol);
density=0.04;
parameters.scale=density/real(z'*z);
parameters.duration=[2;3;5]*1e-3;
parameters.steps=[12;18;30];
parameters.acquire=true(3,1);
parameters.gradient=[0 0 1;0 0 -0.5;0 0 0]*1e-3;
parameters.rf=[300+200i;-100+40i;0];
fid=imaging(s,@ddf_sequence,parameters);

% Integrate the torque of the local dipolar field at the selected voxel
gamma=spin('1H');
coord=(8-(parameters.npts(3)+1)/2)*parameters.dims(3)/parameters.npts(3);
dipolar=gamma*1.25663706127e-6*density;
moment=[0;0;1];
expected=zeros(3,1);
for n=1:3
    external=[real(parameters.rf(n));imag(parameters.rf(n));-gamma*coord*parameters.gradient(n,3)];
    rhs=@(~,m)cross(external+[0;0;dipolar*m(3)],m);
    [~,trajectory]=ode45(rhs,[0 parameters.duration(n)],moment,...
                         odeset('AbsTol',1e-12,'RelTol',1e-10));
    moment=trajectory(end,:).';
    expected(n)=real(z'*z)*(moment(1)-1i*moment(2));
end
result=test_close(result,'three waveform edges',fid,expected,2e-5,0,...
                  'the axial tensor reduces to an Mz-dependent offset for any tilt');

% Isolate resonant uniform-coil radiation damping from the dipolar field
parameters.rho0_st={coil_state(s,'Lx',1)};
parameters.gradient=zeros(1,3);
parameters.rf=0;
parameters.acquire=true;
parameters.duration=0.04;
parameters.steps=160;
damping=1/(gamma*density*parameters.duration);
parameters.field=@(m)repmat(damping*[-mean(m(:,2)) mean(m(:,1)) 0],size(m,1),1);
damped=imaging(s,@ddf_sequence,parameters);
parameters.scale=0;
control=imaging(s,@ddf_sequence,parameters);
result=test_close(result,'uniform-coil radiation damping',damped/control,sech(1),...
                  1e-6,0,'an initially transverse moment relaxes as sech(t/T_rd)');

end


