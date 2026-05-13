% Tests frame-transformation and averaging front-end kernels. Syntax:
%
%                    result=test_dynamic_frame_frontends()
%
% Outputs:
%
%     result  - regression test result with explanatory messages
%
% The test exercises carrier(), frqoffset(), rotframe(), average(), and
% orientation() on compact systems and analytically controlled limiting
% cases.
%
% ilya.kuprov@weizmann.ac.il

function result=test_dynamic_frame_frontends()

% Announce the test target
fprintf('TESTING: Dynamic frame and averaging front ends\n');

% State the dynamic frame target of the test
result=new_test_result('kernel/dynamic_frame_frontends',...
                       'Dynamic frame and averaging front ends',...
                       'frame-transformation helpers must preserve exact algebraic limits.');

% Check carrier Hamiltonian operator-type algebra
result=local_test_carrier(result);

% Check single-channel and duplicate-channel frequency offsets
result=local_test_frqoffset(result);

% Check the zeroth-order rotating-frame transformation
result=local_test_rotframe(result);

% Check average-Hamiltonian and Krylov-Bogolyubov branches
result=local_test_average(result);

% Check explicit non-zero-angle rotational-basis contraction
result=local_test_orientation(result);

end


function result=local_test_carrier(result)

% Build a heteronuclear Liouville-space system with non-zero carriers
spin_system=local_sphten_system();
basefrqs=spin_system.inter.basefrqs;

% Evaluate all Liouville-space operator-type paths for the proton carrier
C_left=carrier(spin_system,'1H','left');
C_right=carrier(spin_system,'1H','right');
C_comm=carrier(spin_system,'1H','comm');
C_acomm=carrier(spin_system,'1H','acomm');

% Compare commutator and anticommutator algebra against side products
result=test_close(result,'carrier commutator algebra',C_comm,C_left-C_right,1e-8,1e-14,...
                  'carrier() comm mode must equal left-minus-right side products');
result=test_close(result,'carrier anticommutator algebra',C_acomm,C_left+C_right,1e-8,1e-14,...
                  'carrier() acomm mode must equal left-plus-right side products');

% Compare isotope-specific and all-spin carrier selection
C_all=carrier(spin_system,'all','comm');
C_ref=carrier(spin_system,'1H','comm')+carrier(spin_system,'13C','comm');
result=test_close(result,'carrier all-spin selection',C_all,C_ref,1e-8,1e-14,...
                  'carrier() all-spin mode must sum the isotope-specific carriers');

% Compare the proton carrier against the corresponding Lz operator
C_ref=basefrqs(1)*operator(spin_system,{'Lz'},{1},'comm');
result=test_close(result,'carrier proton frequency',C_comm,C_ref,1e-8,1e-14,...
                  'carrier() must scale Lz by the free-particle Larmor frequency');

end


function result=local_test_frqoffset(result)

% Build the offset test system and start from a zero Hamiltonian
spin_system=local_sphten_system();
H0=sparse(size(spin_system.bas.basis,1),size(spin_system.bas.basis,1));

% Apply independent offsets to both isotopes
parameters.spins={'1H','13C'};
parameters.offset=[25.0,-40.0];
H_obs=frqoffset(spin_system,H0,parameters);
H_ref=2*pi*25.0*operator(spin_system,'Lz','1H')-...
      2*pi*40.0*operator(spin_system,'Lz','13C');
result=test_close(result,'frequency offset independent channels',H_obs,H_ref,1e-12,1e-12,...
                  'frqoffset() must add one Lz term per independent channel');

% Apply duplicated proton channels carrying the same offset
parameters.spins={'1H','13C','1H'};
parameters.offset=[25.0,-40.0,25.0];
H_obs=frqoffset(spin_system,H0,parameters);
result=test_close(result,'frequency offset duplicate channels',H_obs,H_ref,1e-12,1e-12,...
                  'frqoffset() must merge duplicate channel names with equal offsets');

end


function result=local_test_rotframe(result)

% Build a one-spin Hilbert-space laboratory-frame system
sys.magnet=14.1;
sys.isotopes={'1H'};
inter.zeeman.scalar={0.0};
bas.formalism='zeeman-hilb';
bas.approximation='none';
spin_system=test_spin_system(sys,inter,bas);
spin_system=assume(spin_system,'labframe');

% Form a carrier plus a small transverse perturbation
H0=carrier(spin_system,'1H');
H1=1e3*operator(spin_system,'Lx','1H');
H1=(H1+H1')/2;

% Zeroth-order rotating-frame theory must remove only the carrier
Hr=rotframe(spin_system,H0,H0+H1,'1H',0);
result=test_close(result,'rotframe zeroth-order perturbation',Hr,H1,1e-10,1e-14,...
                  'rotframe() at order zero must return the laboratory-frame perturbation H-H0');

end


function result=local_test_average(result)

% Build a quiet spin system for averaging diagnostics
spin_system=local_sphten_system();

% Define an unmodulated two-level Hamiltonian decomposition
Hp=sparse(2,2);
H0=sparse([0 1;1 0]);
Hm=sparse(2,2);
omega=2*pi*1000;
theories={'ah_first_order','ah_second_order','ah_third_order',...
          'kb_first_order','kb_second_order','kb_third_order','matrix_log'};

% Check every advertised averaging branch in the unmodulated limit
for n=1:numel(theories)
    H_obs=average(spin_system,Hp,H0,Hm,omega,theories{n});
    result=test_close(result,['average unmodulated ' theories{n}],H_obs,H0,1e-10,1e-12,...
                      'average() must leave an unmodulated Hamiltonian unchanged');
end

end


function result=local_test_orientation(result)

% Build a sparse rank-one and rank-two rotational basis by hand
Q=cell(2,1);
Q{1}=cell(3,3);
Q{2}=cell(5,5);
for row=1:3
    for col=1:3
        Q{1}{row,col}=sparse(2,2);
    end
end
for row=1:5
    for col=1:5
        Q{2}{row,col}=sparse(2,2);
    end
end
Q{1}{1,3}=sparse([0 1;0 0]);
Q{1}{3,1}=sparse([0 0;1 0]);
Q{2}{3,3}=sparse([1 0;0 -1]);

% Contract the basis at a non-trivial orientation
angles=[0.2 0.3 -0.4];
H_obs=orientation(Q,angles);
D1=wigner(1,angles(1),angles(2),angles(3));
D2=wigner(2,angles(1),angles(2),angles(3));
H_ref=D1(1,3)*Q{1}{1,3}+D1(3,1)*Q{1}{3,1}+D2(3,3)*Q{2}{3,3};
H_ref=(H_ref+H_ref')/2;

% Check the non-zero-angle sparse contraction
result=test_close(result,'non-zero Euler orientation',H_obs,H_ref,1e-14,1e-14,...
                  'orientation() must contract populated Wigner components at non-zero Euler angles');

end


function spin_system=local_sphten_system()

% Build the common heteronuclear spherical-tensor test system
sys.magnet=14.1;
sys.isotopes={'1H','13C'};
inter.zeeman.scalar={1.0,2.0};
inter.coupling.scalar=cell(2);
inter.coupling.scalar{1,2}=10.0;
inter.coupling.scalar{2,2}=0.0;
bas.formalism='sphten-liouv';
bas.approximation='none';
spin_system=test_spin_system(sys,inter,bas);

end


