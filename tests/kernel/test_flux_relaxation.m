% Explicit flux relaxation regression. Syntax:
%
%                    result=test_flux_relaxation()
%
% Outputs:
%
%    result - flux-switch regression assertions and failures
%
% ilya.kuprov@weizmann.ac.il

function result=test_flux_relaxation()

% Configure a small directed-transfer system
result=new_test_result('kernel/flux_relaxation','Flux relaxation switch',...
                       'Activation, additivity, ownership and validation.');
sys.magnet=9.4; sys.isotopes={'1H','1H'}; sys.parallel={'local',1};
sys.output='hush'; sys.disable={'hygiene'}; sys.parprops={};
inter.zeeman.scalar={0,0}; inter.flux_rate=[0 3;7 0];
inter.flux_type='intramolecular'; inter.rlx_keep='labframe';
inter.equilibrium='zero'; inter.relaxation={'flux'};
bas.formalism='sphten-liouv'; bas.approximation={'none'};
loud=sys; loud.output='console';
summary=evalc('s=create(loud,inter);'); s=basis(s,bas);
F=rlx_flux(s); s.sys.output='hush';

% Check a directed column against the two-site transfer equation
source=state(s,'Lz',1); destination=state(s,'Lz',2);
expected=3*(destination-source);
result=test_close(result,'independent directed transfer',F*source,expected,...
                  100*eps(norm(expected,1)),0,...
                  'A unit source loses rate three and the other spin gains it.');

result=test_true(result,'flux active',nnz(F)>0,'Flux must transfer magnetisation.');
result=test_close(result,'flux alone',relaxation(s),F,100*eps(norm(F,1)),0,...
                  'Flux alone equals the standalone generator.');
result=test_true(result,'internal ownership',...
                 isequal(s.rlx.flux_rate,inter.flux_rate)&&...
                 strcmp(s.rlx.flux_type,inter.flux_type),...
                 'Inputs compile into rlx fields.');
result=test_true(result,'enabled summary',contains(lower(summary),'flux'),...
                 'Enabled flux must be reported.');

% Stored rates do not enable the mechanism
without=rmfield(inter,'relaxation');
quiet=evalc('off=create(loud,without);'); off=basis(off,bas);
result=test_close(result,'flux inactive without switch',relaxation(off),...
                  sparse(size(F,1),size(F,2)),0,0,...
                  'Flux parameters alone are inert.');
result=test_close(result,'kinetics excludes flux',kinetics(s),kinetics(off),...
                  0,0,'Flux must not enter kinetics.');
result=test_true(result,'disabled summary absent',...
                 ~contains(lower(quiet),'magnetisation flux'),...
                 'Disabled flux must not print a summary.');

% The same inactive inputs remain safe in the Zeeman formalism
zeeman=bas; zeeman.formalism='zeeman-liouv';
off_z=test_spin_system(sys,without,zeeman);
result=test_true(result,'disabled flux in Zeeman',nnz(relaxation(off_z))==0,...
                 'Inactive nonzero flux metadata does not invoke its formalism gate.');

% Combine flux with two independent relaxation mechanisms
t1=without; t1.relaxation={'t1_t2'};
t1.r1_rates={0.5,1/3}; t1.r2_rates={1,2/3};
t=create(sys,t1); t=basis(t,bas); T=relaxation(t);
t1.relaxation={'t1_t2','flux'};
tf=create(sys,t1); tf=basis(tf,bas);
result=test_close(result,'T1 T2 flux additive',relaxation(tf),T+F,...
                  100*eps(norm(F,1)+norm(T,1)),0,'Flux contribution appears exactly once.');

% Include anisotropic interactions so the Redfield baseline is nonzero
red=without; red.relaxation={'redfield'}; red.tau_c={1e-9};
red.zeeman.matrix={diag([10 15 20]),diag([5 10 15])};
r=create(sys,red); r=basis(r,bas); R=relaxation(r);
result=test_true(result,'nonzero Redfield baseline',nnz(R)>0,...
                 'Combined-mechanism additivity must exercise actual relaxation.');
red.relaxation={'redfield','flux'};
rf=create(sys,red); rf=basis(rf,bas);
result=test_close(result,'Redfield flux additive',relaxation(rf),R+F,...
                  100*eps(norm(F,1)+norm(R,1)),0,'Flux is included once alongside Redfield.');
zero=inter; zero.flux_rate=zeros(2);
z=create(sys,zero); z=basis(z,bas);
result=test_close(result,'zero flux',relaxation(z),...
                  sparse(size(F,1),size(F,2)),0,0,...
                  'Zero flux rates must be valid.');
zero.relaxation={'t1_t2','flux'}; zero.r1_rates={0.5,1/3}; zero.r2_rates={1,2/3};
z=create(sys,zero); z=basis(z,bas);
result=test_close(result,'zero flux with T1 T2',relaxation(z),T,...
                  100*eps(norm(T,1)),0,'Zero flux preserves other models.');

% Reject incomplete new inputs and the retired chemistry-owned syntax
invalid={rmfield(inter,'flux_rate'),rmfield(inter,'flux_type')};
legacy=rmfield(inter,{'flux_rate','flux_type'});
legacy.relaxation={};
legacy.chem.flux_rate=inter.flux_rate;
legacy.chem.flux_type=inter.flux_type;
invalid{end+1}=legacy;
for n=1:numel(invalid)
    rejected=false;
    try
        create(sys,invalid{n});
    catch err
        if n==3
            rejected=strcmp(err.identifier,'Spinach:create:chemistryField');
        else
            rejected=contains(err.message,'inter.flux_');
        end
    end
    result=test_true(result,['invalid inputs ' num2str(n)],rejected,...
                     'Missing paired or obsolete inputs must fail.');
end
part=sys; % Merge and spin deletion retain the top-level input and rlx storage contract
part.isotopes={'1H'};
a=inter; a.flux_rate=0; a.zeeman.scalar={0};
[~,merged]=merge_inp({part,part},{a,a});
result=test_true(result,'merge top-level fields',...
                 isfield(merged,'flux_rate')&&isfield(merged,'flux_type')&&...
                 isequal(merged.flux_rate,zeros(2)),...
                 'Merge keeps top-level flux configuration.');
b=a; b.flux_rate=[0 1;2 0]; rejected=false;
try, merge_inp({part,part},{a,b}); catch, rejected=true; end
result=test_true(result,'merge shape validation',rejected,...
                 'Subsystem rate shape must match its spins.');
reduced=kill_spin(s,2);
result=test_true(result,'kill spin reindex',...
                 isequal(reduced.rlx.flux_rate,0),...
                 'Deleted spin removes rate row and column.');
end

