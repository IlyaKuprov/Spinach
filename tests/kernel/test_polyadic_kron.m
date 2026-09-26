% Tests affixed polyadic tensor products and spin-space flow lifting. Syntax:
%
%                    result=test_polyadic_kron()
%
% Outputs:
%
%     result - regression results against explicit matrix references
%
% Complex rectangular controls cover both operand orders and multiple
% affixes. Production v2fplanck calls compare scalar, voxel, and shear
% velocities in sparse and polyadic representations.
%
% talos@spindynamics.org

function result=test_polyadic_kron()

% State the target and preserve the caller's random stream
result=new_test_result('kernel/polyadic_kron','Affixed polyadic tensor products',...
                       'Kronecker products must preserve all prefixes and suffixes.');
old_rng=rng();
restore_rng=onCleanup(@()rng(old_rng));
rng(240924);

% Build complex non-Hermitian sums with rectangular multiple affixes
for trial=1:3
    a=randn(2)+1i*randn(2);
    b=randn(3)+1i*randn(3);
    c=randn(2)+1i*randn(2);
    d=randn(3)+1i*randn(3);
    ext=randn(2,3)+1i*randn(2,3);
    left1=sparse(randn(5,6)+1i*randn(5,6));
    left2=sparse(randn(7,5)+1i*randn(7,5));
    right1=sparse(randn(6,4)+1i*randn(6,4));
    right2=sparse(randn(4,8)+1i*randn(4,8));
    P=polyadic({{a,b},{c,d}});
    base=kron(a,b)+kron(c,d);
    objects={P,left1*P,left2*(left1*P),P*right1,P*right1*right2,...
             left2*(left1*P)*right1*right2};
    refs={base,left1*base,left2*left1*base,base*right1,...
          base*right1*right2,left2*left1*base*right1*right2};
    assert(numel(objects{3}.prefix)==2&&numel(objects{5}.suffix)==2);

    % Compare full expansion and complex two-column action in both orders
    for shape=1:numel(objects)
        for order=1:2
            if order==1
                Q=kron(objects{shape},ext);
                ref=kron(refs{shape},ext);
            else
                Q=kron(ext,objects{shape});
                ref=kron(ext,refs{shape});
            end
            rhs=randn(size(ref,2),2)+1i*randn(size(ref,2),2);
            label=sprintf('trial %d shape %d order %d',trial,shape,order);
            result=test_close(result,[label ' full'],full(Q),ref,1e-12,1e-12,...
                              'Tensor extension must preserve every matrix factor.');
            result=test_close(result,[label ' action'],Q*rhs,ref*rhs,1e-12,1e-12,...
                              'Deferred action must agree with explicit expansion.');
        end
    end
end

% Build a physical one-spin Liouville system with one process worker
sys.magnet=14.1;
sys.isotopes={'1H'};
sys.parallel={'processes',1};
sys.enable={'polyadic'};
sys.disable={'hygiene'};
sys.output='hush';
inter.zeeman.scalar={0};
spin_system=create(sys,inter);
bas.formalism='sphten-liouv';
bas.approximation='none';
spin_system=basis(spin_system,bas);
dense_system=spin_system;
dense_system.sys.enable={};

% Compare scalar flow with identical voxel flow and divergence-free shear
parameters.npts=[10 10];
parameters.dims=[0.01 0.01];
parameters.deriv={'period',3};
parameters.v=0;
velocities={1e-3,1e-3*ones(10),ones(10,1)*(1:10)*1e-3};
for flow_case=1:numel(velocities)
    parameters.u=velocities{flow_case};
    ref=full(v2fplanck(dense_system,parameters));
    Q=v2fplanck(spin_system,parameters);
    if flow_case==1
        uniform_ref=ref;
    elseif flow_case==2
        result=test_close(result,'scalar versus voxel flow',ref,uniform_ref,1e-12,1e-12,...
                          'Constant voxel velocities describe the same uniform flow.');
    end
    if ~isscalar(parameters.u)
        [flow_x,~,~]=hydrodynamics(spin_system,parameters);
        result=test_close(result,'flow divergence',flow_x*parameters.u(:),zeros(100,1),...
                          1e-12,1e-12,'The chosen transverse shear is divergence-free.');
    end
    rhs=randn(size(ref,2),2)+1i*randn(size(ref,2),2);
    label=sprintf('flow %d',flow_case);
    result=test_close(result,[label ' dimensions'],double(size(Q)),[400 400],0,0,...
                      'Spatial flow must lift into the four-dimensional spin space.');
    result=test_close(result,[label ' full'],full(Q),ref,1e-12,1e-12,...
                      'Polyadic and sparse production generators must agree.');
    result=test_close(result,[label ' action'],Q*rhs,ref*rhs,1e-12,1e-12,...
                      'Spin-space flow action must not depend on its representation.');
end

end


