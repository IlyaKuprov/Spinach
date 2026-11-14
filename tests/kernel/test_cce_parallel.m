% Deterministic serial/parallel scheduling and resource acceptance test.
%
% ilya.kuprov@weizmann.ac.il

function result=test_cce_parallel()

% Declare the acceptance test and report unavailable parallel execution
result=new_test_result('kernel/cce_parallel','CCE deterministic scheduling',...
                       'Serial, one/two workers, common samples and resource caps.');
if ~license('test','Distrib_Computing_Toolbox')
    error('test_cce_parallel:unavailable',...
          'SKIPPED: Parallel Computing Toolbox unavailable; this harness has no skip status.');
end

% Construct a small analytic Ising bath with deterministic sampled states
sys.isotopes={'13C','13C','13C'}; sys.magnet=0; sys.output='hush';
sys.disable={'hygiene'}; sys.parallel={'local',1};
inter.zeeman.matrix={zeros(3),zeros(3),zeros(3)};
inter.coupling.matrix=cell(3); inter.coupling.matrix{1,2}=diag([0 0 17]);
inter.coupling.matrix{1,3}=diag([0 0 29]);
spin_system=create(sys,inter);
p=struct('central_spins',1,'bath_spins',[2 3],'method','cce',...
         'expansion','standard','quantity','coherence','times',[0 0.001 0.004],...
         'order',2,'central',struct('alpha',[1;0],'beta',[0;1]),...
         'assumptions','labframe','orientation',[0 0 0],'interaction_cutoff',Inf,...
         'graph',struct('mode','complete','score','coupling_norm'),...
         'meanfield','none','sequence',[],'averaging','after_expansion',...
         'bath_state',struct('mode','sampled','basis','uniform','probabilities','uniform',...
                             'samples',4,'seed',193),...
         'numerics',struct('atol',1e-8,'rtol',1e-8,'max_memory_bytes',2^26));

% Measure serial and one-worker execution on the identical sample list
p.parallel=struct('mode','serial','batch_size',3,'max_tasks',16);
started=tic; serial=cce(spin_system,p); serial_time=toc(started);
p.parallel.mode='clusters'; started=tic;
one=cce(spin_system,p); one_time=toc(started);
result=test_true(result,'one-worker raw identity',isequaln(serial.raw,one.raw),...
                 'Indexed raw records are bitwise identical.');
result=test_close(result,'one-worker reduction',one.value,serial.value,0,0,...
                  'The final reduction runs only on the client.');

% Request two workers through create rather than a hidden scheduler pool
sys.parallel={'local',2}; spin_system=create(sys,inter);
started=tic; two=cce(spin_system,p); two_time=toc(started);
result=test_true(result,'two-worker sample identity',isequal(serial.ensemble.assignments,two.ensemble.assignments),...
                 'State substreams depend on seed and state ID, never worker rank.');
result=test_true(result,'two-worker raw identity',isequaln(serial.raw,two.raw),...
                 'Cluster records have a deterministic gather order.');
result=test_close(result,'two-worker reduction',two.value,serial.value,0,0,...
                  'The complex result is bitwise identical across worker counts.');
p.parallel.mode='states'; states=cce(spin_system,p);
result=test_close(result,'state scheduling',states.value,serial.value,0,0,...
                  'There is one outer parallel layer for either scheduling mode.');

% Compare against independent analytic conditional Ising phases
assignments=serial.ensemble.assignments; magnetic=1.5-assignments;
reference=mean(exp(-1i*(magnetic*[17;29])*2*pi.*p.times),1).';
result=test_close(result,'analytic sampled coherence',serial.value,reference,1e-12,1e-12,...
                  'Conditional phases follow the exact diagonal Ising Hamiltonian.');

% Refuse oversized task ledgers and per-worker storage before propagation
p.parallel.max_tasks=15; refused=false;
try
    cce(spin_system,p); catch exception, refused=strcmp(exception.identifier,'cce_schedule:tasks'); end
result=test_true(result,'task cap',refused,'Sixteen tasks cannot fit a fifteen-task cap.');
p.parallel.max_tasks=16; p.numerics.max_memory_bytes=2000; refused=false;
try
    cce(spin_system,p); catch exception, refused=strcmp(exception.identifier,'cce_schedule:memory'); end
result=test_true(result,'memory cap',refused,'Raw storage and dense working matrices are bounded together.');

% Print the pilot timings without imposing an unjustified speedup target
fprintf('CCE pilot serial %.6g s; one worker %.6g s; two workers %.6g s; discrepancy %.3g.\n',...
        serial_time,one_time,two_time,max(abs(serial.value-two.value),[],'all'));
sys.parallel={'local',1}; create(sys,inter);

end


