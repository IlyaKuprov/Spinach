% Tests explicit assumptions throughout rotor-stack construction. Syntax:
%
%                    result=test_rotor_assume()
%
% Outputs:
%
%     result  - regression test result with explanatory messages
%
% Fresh, stale, and matching objects must give identical numerical-frame
% stacks. Both MAS frames and Hilbert/Liouville representations are tested.
%
% talos@spindynamics.org

function result=test_rotor_assume()

% Describe the regression target
result=new_test_result('kernel/rotor_assume','Rotor-stack assumptions',...
                      'Explicit assumptions must reach numerical rotating frames.');

% Specify an anisotropic heteronuclear pair with noncommuting lab Hamiltonians
sys.magnet=9.4;
sys.isotopes={'1H','13C'};
sys.output='hush';
sys.disable={'hygiene'};
sys.parallel={'processes',1};
inter.zeeman.eigs={[-12 5 20],[-30 10 45]};
inter.zeeman.euler={[0.2 0.5 0.7],[0.6 0.4 0.3]};
inter.coupling.scalar=cell(2);
inter.coupling.scalar{1,2}=150;
bas.approximation='none';

% Specify a tilted rotor and a nonzero transmitter offset
parameters.axis=[1 2 3]/sqrt(14);
parameters.offset=[70 -35];
parameters.spins={'1H','13C'};
parameters.max_rank=1;
parameters.orientation=[0.3 0.7 0.2];
formalisms={'zeeman-hilb','sphten-liouv'};
masframes={'rotor','magnet'};

% Exercise both formalisms and both MAS orientation conventions
for f=1:numel(formalisms)
    bas.formalism=formalisms{f};
    spin_system=basis(create(sys,inter),bas);
    matched=assume(spin_system,'labframe');
    stale=assume(spin_system,'nmr');

    % Scale comparison accuracy to the laboratory carriers used by dirdiff
    carrier_norm=norm(carrier(matched,'1H'),'fro')+...
                 norm(carrier(matched,'13C'),'fro');
    frame_tol=1e-14*carrier_norm*sqrt(2*parameters.max_rank+1);
    for m=1:numel(masframes)
        parameters.masframe=masframes{m};
        parameters.rframes={{'1H',1},{'13C',1}};

        % Compare fresh and stale objects against matching explicit assumptions
        [reference,phases]=rotor_stack(matched,parameters,'labframe');
        [fresh,fresh_phases]=rotor_stack(spin_system,parameters,'labframe');
        [old,old_phases]=rotor_stack(stale,parameters,'labframe');
        repeated=rotor_stack(matched,parameters,'labframe');
        result=test_close(result,'repeated matching frames',cat(3,repeated{:}),...
                          cat(3,reference{:}),frame_tol,1e-12,...
                          'repeated matching objects bound numerical-frame repeatability');
        result=test_close(result,'fresh numerical frames',cat(3,fresh{:}),...
                          cat(3,reference{:}),frame_tol,1e-12,...
                          'fresh objects must use the requested laboratory assumptions');
        result=test_close(result,'stale numerical frames',cat(3,old{:}),...
                          cat(3,reference{:}),frame_tol,1e-12,...
                          'previous NMR assumptions must not leak into frame validation');
        result=test_close(result,'fresh rotor phases',fresh_phases,phases,0,0,...
                          'assumption history must not change the rotor phase grid');
        result=test_close(result,'stale rotor phases',old_phases,phases,0,0,...
                          'assumption history must not change the rotor phase grid');

        % Keep rejection of explicitly inconsistent numerical-frame requests
        rejected=false;
        try
            rotor_stack(matched,parameters,'nmr');
        catch failure
            rejected=contains(failure.message,'already in the rotating frame');
        end
        result=test_true(result,'explicit NMR frame rejection',rejected,...
                         'numerical frames require laboratory-frame spins');

        % Check ordinary empty-frame NMR stacks against matching NMR assumptions
        parameters.rframes={};
        ordinary=rotor_stack(spin_system,parameters,'nmr');
        normal_ref=rotor_stack(stale,parameters,'nmr');
        result=test_close(result,'ordinary NMR stack',cat(3,ordinary{:}),...
                          cat(3,normal_ref{:}),0,0,...
                          'ordinary rotor stacks must remain independent of prior assumptions');

        % Verify that the laboratory control is complex and genuinely noncommuting
        lab=rotor_stack(spin_system,parameters,'labframe');
        result=test_true(result,'complex laboratory stack',norm(imag(lab{1}),'fro')>1,...
                         'tilted anisotropic interactions must exercise complex matrices');
        result=test_true(result,'noncommuting laboratory stack',...
                         norm(lab{1}*lab{2}-lab{2}*lab{1},'fro')>1,...
                         'different rotor phases must not share a diagonal eigenbasis');
    end
end

end


