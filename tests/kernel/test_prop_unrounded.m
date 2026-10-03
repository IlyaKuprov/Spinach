% Tests numerical termination of unrounded Taylor propagation. Syntax:
%
%                    result=test_prop_unrounded()
%
% Outputs:
%
%     result - sparse/dense, disabled-cleanup, zero-tolerance,
%              nonnormal, nilpotent, and cached propagation checks
%
% talos@spindynamics.org

function result=test_prop_unrounded()

% State the unrounded-series target
result=new_test_result('kernel/prop_unrounded','Unrounded Taylor termination',...
                       'explicitly stored sparse zeros must not keep an exhausted Taylor series running.');

% Create a quiet system with the default Taylor-branch boundary
sys.magnet=1; sys.isotopes={'1H'}; sys.output='hush';
sys.enable={}; sys.disable={'hygiene'};
sys.parallel={'processes',1}; sys.parprops={};
inter.zeeman.scalar={0}; spin_system=create(sys,inter);
spin_system.tols.small_matrix=200; dt=0.13;
L=spdiags([-0.31*ones(200,1) (1:200)'/200 0.73*ones(200,1)],-1:1,200,200);
L=L-0.2i*speye(200); reference=expm(full(-1i*L*dt));

% Test both ways of disabling term chopping with sparse and dense input
for policy=1:2
    spin_system.sys.disable={'hygiene'}; spin_system.tols.prop_chop=1e-10;
    if policy==1
        spin_system.sys.disable={'hygiene','clean-up'};
    else
        spin_system.tols.prop_chop=0;
    end
    for representation=1:2
        generator=L;
        if representation==2, generator=full(generator); end
        P=propagator(spin_system,generator,dt);
        result=test_close(result,'unrounded exponential',P,reference,1e-12,1e-12,...
                          'nonnormal dissipative propagation must terminate and agree with expm');
    end
end

% Check zero, negative, and imaginary timesteps without chopping
spin_system.tols.prop_chop=0;
for step=[0 -0.13 0.13i 3.1]
    P=propagator(spin_system,L,step);
    result=test_close(result,'complex or signed timestep',P,expm(full(-1i*L*step)),1e-12,1e-12,...
                      'termination must preserve the documented finite complex timestep contract');
end

% Retain the chopped normal path and a finitely terminating nilpotent control
spin_system.tols.prop_chop=1e-10;
P=propagator(spin_system,L,dt);
result=test_close(result,'chopped nonnormal exponential',P,reference,1e-8,1e-8,...
                  'the normal chopped series must retain its configured accuracy');
spin_system.sys.disable={'hygiene','clean-up'};
L=spdiags((1:200)'/200,1,200,200);
P=propagator(spin_system,L,dt);
result=test_close(result,'nilpotent exponential',P,expm(full(-1i*L*dt)),1e-12,1e-12,...
                  'a finite nonnormal Taylor series must retain its exact stopping behaviour');

% Exercise a fresh shared-cache miss followed by a repeated hit
spin_system.sys.enable={'prop_cache'};
L=spdiags([-0.31*ones(200,1) (1:200)'/200 0.73*ones(200,1)],-1:1,200,200);
L=L-0.2i*speye(200); P=propagator(spin_system,L,dt); Q=propagator(spin_system,L,dt);
result=test_close(result,'unrounded cached miss',P,reference,1e-12,1e-12,...
                  'a cached unrounded miss must terminate before its result can be inserted');
result=test_close(result,'unrounded cached hit',Q,P,0,0,...
                  'a repeated cache hit must reproduce the terminating miss');

end


