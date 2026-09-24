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
% Already-rotating spins are rejected across all frame assumption sets,
% while laboratory nuclei in mixed-frame systems remain transformable.
% Carrier-free solid-effect components retain valid empty-frame stacks,
% but direct and fresh/stale numerical-frame requests are refused.
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

% Specify a hyperfine-coupled pair with laboratory-frame nuclear dynamics
sys.magnet=0.35;
sys.isotopes={'E','1H'};
inter.zeeman.eigs={[2.0023 2.0023 2.0023],[-12 5 20]};
inter.coupling.scalar=cell(2);
inter.coupling.eigs=cell(2);
inter.coupling.euler=cell(2);
inter.coupling.eigs{1,2}=[1e4 2e4 4e4];
inter.coupling.euler{1,2}=[0.4 0.6 0.8];
parameters.spins={'E','1H'};
parameters.masframe='rotor';
assumptions={'nmr','cavity','esr','deer','deer-zz','spin-phonon'};

% Cover every all-spin and electron-only rotating-frame assumption set
for f=1:numel(formalisms)
    bas.formalism=formalisms{f};
    spin_system=basis(create(sys,inter),bas);
    stale=assume(spin_system,'labframe');
    for a=1:numel(assumptions)
        assumed=assume(spin_system,assumptions{a});
        parameters.rframes={};
        ordinary=rotor_stack(spin_system,parameters,assumptions{a});
        normal_ref=rotor_stack(assumed,parameters,assumptions{a});
        result=test_close(result,[assumptions{a} ' ordinary stack'],...
                          cat(3,ordinary{:}),cat(3,normal_ref{:}),0,0,...
                          'empty frame requests must remain valid');
        targets={'E'};
        if ismember(assumptions{a},{'nmr','cavity'})
            targets={'E','1H'};
        end
        for t=1:numel(targets)

            % Reject direct numerical transformations of already-rotating spins
            H0=carrier(assumed,targets{t});
            H=(ordinary{1}+ordinary{1}')/2;
            rejected=false;
            try
                rotframe(assumed,H0,H,targets{t},1);
            catch failure
                rejected=contains(failure.message,'already in the rotating frame');
            end
            result=test_true(result,[assumptions{a} ' direct ' targets{t}],rejected,...
                             'rotframe must reject a second carrier subtraction');

            % Reject the same request from fresh and stale rotor-stack inputs
            parameters.rframes={{targets{t},1}};
            inputs={spin_system,stale};
            for s=1:numel(inputs)
                rejected=false;
                try
                    rotor_stack(inputs{s},parameters,assumptions{a});
                catch failure
                    rejected=contains(failure.message,'already in the rotating frame');
                end
                result=test_true(result,[assumptions{a} ' rotor ' targets{t}],rejected,...
                                 'explicit assumptions must reject fresh and stale inputs alike');
            end
        end
    end

    % Retain nuclear numerical frames under all electron-only rotating sets
    parameters.rframes={{'1H',1}};
    reference=rotor_stack(spin_system,parameters,'esr');
    for a=3:numel(assumptions)
        fresh=rotor_stack(spin_system,parameters,assumptions{a});
        old=rotor_stack(stale,parameters,assumptions{a});
        frame_tol=1e-14*norm(carrier(spin_system,'1H'),'fro')*sqrt(numel(fresh));
        result=test_close(result,[assumptions{a} ' nuclear frame'],...
                          cat(3,fresh{:}),cat(3,reference{:}),frame_tol,1e-12,...
                          'one-electron systems share the ESR nuclear-frame Hamiltonian');
        result=test_close(result,[assumptions{a} ' stale nuclear frame'],...
                          cat(3,old{:}),cat(3,reference{:}),frame_tol,1e-12,...
                          'laboratory nuclear spins must remain eligible for transformation');
    end

    % Verify that the mixed-frame control is complex and noncommuting
    parameters.rframes={};
    mixed=rotor_stack(spin_system,parameters,'spin-phonon');
    result=test_true(result,'complex mixed-frame stack',norm(imag(mixed{1}),'fro')>1,...
                     'pseudosecular hyperfine terms must exercise complex matrices');
    result=test_true(result,'noncommuting mixed-frame stack',...
                     norm(mixed{1}*mixed{2}-mixed{2}*mixed{1},'fro')>1,...
                     'rotor phases must retain noncommuting nuclear dynamics');
end

% Refuse numerical frames on carrier-free solid-effect Hamiltonian components
inter.coupling.eigs{1,2}=[0 0 0];
parameters.offset=[0 0];
components={'se_dnp_h+','se_dnp_h-','se_dnp_h0'};
for f=1:numel(formalisms)
    bas.formalism=formalisms{f};
    spin_system=basis(create(sys,inter),bas);
    stale=assume(spin_system,'labframe');
    for a=1:numel(components)
        assumed=assume(spin_system,components{a});
        parameters.rframes={};
        ordinary=rotor_stack(spin_system,parameters,components{a});
        result=test_true(result,[components{a} ' zero component'],...
                         all(cellfun(@(H) nnz(H)==0,ordinary)),...
                         'uncoupled solid-effect components contain no Zeeman carrier');
        for isotope={'E','1H'}
            H0=carrier(assumed,isotope{1});
            rejected=false;
            try
                rotframe(assumed,H0,ordinary{1},isotope{1},1);
            catch failure
                rejected=contains(failure.message,'solid-effect Hamiltonian components');
            end
            result=test_true(result,[components{a} ' direct ' isotope{1}],rejected,...
                             'a carrier-free component is not a laboratory Hamiltonian');
            parameters.rframes={{isotope{1},1}};
            for input={spin_system,stale}
                rejected=false;
                try
                    rotor_stack(input{1},parameters,components{a});
                catch failure
                    rejected=contains(failure.message,'solid-effect Hamiltonian components');
                end
                result=test_true(result,[components{a} ' rotor ' isotope{1}],rejected,...
                                 'fresh and stale inputs must not acquire a false carrier');
            end
        end
    end
end

% Specify a quadrupolar nucleus alongside an already-rotating spin-half nucleus
parameters.offset=[70 -35];
sys.magnet=9.4;
sys.isotopes={'1H','14N'};
inter.zeeman.eigs={[-12 5 20],[-30 10 45]};
inter.coupling.eigs=cell(2);
inter.coupling.euler=cell(2);
inter.coupling.eigs{2,2}=[-1e5 -2e5 3e5];
inter.coupling.euler{2,2}=[0.4 0.6 0.8];
parameters.spins={'1H','14N'};
for f=1:numel(formalisms)
    bas.formalism=formalisms{f};
    spin_system=basis(create(sys,inter),bas);
    assumed=assume(spin_system,'qnmr');
    stale=assume(spin_system,'nmr');
    parameters.rframes={{'14N',1}};
    reference=rotor_stack(assumed,parameters,'qnmr');
    fresh=rotor_stack(spin_system,parameters,'qnmr');
    old=rotor_stack(stale,parameters,'qnmr');
    frame_tol=1e-14*norm(carrier(assumed,'14N'),'fro')*sqrt(numel(fresh));
    result=test_close(result,'quadrupolar laboratory frame',cat(3,fresh{:}),...
                      cat(3,reference{:}),frame_tol,1e-12,...
                      'qnmr must permit numerical frames for spin-one nuclei');
    result=test_close(result,'stale quadrupolar laboratory frame',cat(3,old{:}),...
                      cat(3,reference{:}),frame_tol,1e-12,...
                      'stale NMR assumptions must not reject a laboratory-frame nucleus');

    % Preserve refusal of a spin-half numerical frame under qnmr
    parameters.rframes={{'1H',1}};
    rejected=false;
    try
        rotor_stack(spin_system,parameters,'qnmr');
    catch failure
        rejected=contains(failure.message,'already in the rotating frame');
    end
    result=test_true(result,'qnmr spin-half rejection',rejected,...
                     'qnmr already removes the spin-half nuclear carrier');
end

end


