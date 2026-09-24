% Tests retained mode data after removing spins or bosonic modes. Syntax:
%
%                    result=test_mode_removal()
%
% Outputs:
%
%     result  - regression checks against independently created systems
%
% The cases include pair couplings, nested first and second derivatives,
% spin-one quadratic terms, complex transverse operators, and mode decay.
%
% talos@spindynamics.org

function result=test_mode_removal()

% State the particle-removal target
result=new_test_result('kernel/mode_removal','Retained mode reindexing',...
                       'particle removal must preserve retained mode Hamiltonians and decay.');

% Cover spectator, spin, mode, simultaneous, logical, and no-op removals
hit_lists={2,3,1,[2 4],[false true false false false],5,[]};
formalisms={'zeeman-hilb','zeeman-liouv'};
for n=1:numel(hit_lists)

    % Reconstruct the retained physical system from its particle identities
    keep=1:5; keep(hit_lists{n})=[];
    original=build_system(1:5);
    trimmed=kill_spin(original,hit_lists{n});
    reference=build_system(keep);
    label=['case ' num2str(n)];
    result=test_true(result,[label ' mode data'],...
                     isequal(trimmed.inter.modes,reference.inter.modes),...
                     'all particle indices, including nested spin leaves, must match');

    % Rebuild the basis and assumptions required by kill_spin
    for k=1:numel(formalisms)
        bas.formalism=formalisms{k}; bas.approximation='none';
        observed=assume(basis(trimmed,bas),'labframe');
        expected=assume(basis(reference,bas),'labframe');
        [H_obs,Q_obs]=hamiltonian(observed);
        [H_ref,Q_ref]=hamiltonian(expected);
        H_obs=H_obs+orientation(Q_obs,[0.17 0.31 0.23]);
        H_ref=H_ref+orientation(Q_ref,[0.17 0.31 0.23]);
        result=test_close(result,[label ' ' formalisms{k}],H_obs,H_ref,...
                          1e-10,1e-12,'retained interactions must match independent reconstruction');

        % Prove that transverse and noncommuting contributions are exercised
        if (n==1)&&(k==1)
            result=test_true(result,'complex Hamiltonian',norm(imag(H_ref),'fro')>1,...
                             'transverse y fields must produce genuinely complex matrix elements');
            mode_num=operator(expected,'N',1);
            result=test_true(result,'noncommuting mode terms',...
                             norm(H_ref*mode_num-mode_num*H_ref,'fro')>1,...
                             'mode couplings and modulation must not commute with mode energy');
            result=test_close(result,'Hermitian Hamiltonian',H_ref,H_ref',1e-10,1e-12,...
                              'real tensor and field derivatives must give a Hermitian Hamiltonian');
        end

        % Compare finite-temperature dissipators in Liouville space
        if k==2
            R_obs=relaxation(observed); R_ref=relaxation(expected);
            result=test_true(result,[label ' nonzero decay'],norm(R_ref,'fro')>1,...
                             'the decay comparison must have a nonzero reference');
            result=test_close(result,[label ' mode decay'],R_obs,R_ref,1e-10,1e-12,...
                              'retained damping and dephasing must match independent reconstruction');
        end
    end
end

end

% Builds the specified physical particles without trimming any input arrays
function spin_system=build_system(keep)

% Check consistency
grumble(keep);

% Assign particle identities and isolated execution settings
isotopes={'C3','13C','2H','1H','V3'};
sys.isotopes=isotopes(keep); sys.magnet=0;
sys.output='hush'; sys.disable={'hygiene'};
sys.parallel={'processes',1}; sys.parprops={};
inter.temperature=1e-7;
count=numel(keep);

% Set distinct mode frequencies, carriers, anharmonicities, and decay rates
fields={'frqs','carriers','anharms','linewidths','t2_times'};
values=[800 1300;10000 17000;13 29;11 19;0.003 0.007];
for n=1:numel(fields)
    inter.modes.(fields{n})=cell(1,count);
    for k=1:2
        index=find(keep==1+4*(k-1));
        if ~isempty(index)
            inter.modes.(fields{n}){index}=values(n,k);
        end
    end
end

% Insert the four pair channels by physical particle identity
fields={'exchange','kerr','longitudinal','dispersive'};
pairs=[1 5;1 5;5 3;1 4]; strengths=[17 19 23 29];
for n=1:numel(fields)
    inter.modes.(fields{n})=cell(count);
    row=find(keep==pairs(n,1)); col=find(keep==pairs(n,2));
    if (~isempty(row))&&(~isempty(col))
        inter.modes.(fields{n}){row,col}=strengths(n);
    end
end

% Build spin-indexed leaves from retained physical spin identities
spin_one=find(keep==3); spin_half=find(keep==4);
tensors=cell(count); vectors=cell(1,count);
if ~isempty(spin_one)
    tensors{spin_one,spin_one}=[7 2 3;2 -5 4;3 4 -2];
    vectors{spin_one}=[13 17 19];
end
if ~isempty(spin_half)
    vectors{spin_half}=[-7 11 3];
end
if (~isempty(spin_one))&&(~isempty(spin_half))
    tensors{spin_one,spin_half}=[3 5 7;11 13 17;19 23 29];
end

% Declare first derivatives, diagonal Raman terms, and mixed Raman terms
inter.modes.coupling_mod=cell(count); inter.modes.zeeman_mod=cell(count);
mode_one=find(keep==1); mode_two=find(keep==5);
for mode=[mode_one mode_two]
    inter.modes.coupling_mod{mode,mode}={tensors,tensors};
    inter.modes.zeeman_mod{mode,mode}={vectors,vectors};
end
if (~isempty(mode_one))&&(~isempty(mode_two))
    inter.modes.coupling_mod{mode_one,mode_two}={[],tensors};
    inter.modes.zeeman_mod{mode_one,mode_two}={[],vectors};
end

% Create a fresh system rather than using the production removal path
spin_system=create(sys,inter);

end

% Consistency enforcement
function grumble(keep)
if (~isnumeric(keep))||(~isrow(keep))||isempty(keep)||...
   any(~ismember(keep,1:5))||any(diff(keep)<=0)
    error('keep must be an increasing row vector of particle identities from 1 to 5.');
end
end


