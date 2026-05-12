% Tests deterministic metadata, hashing, and partition helpers. Syntax:
%
%                    result=test_dynamic_metadata_partition_suite()
%
% Outputs:
%
%     result  - regression test result with explanatory messages
%
% The test checks hashing stability, duplicate-row removal, parallel-state
% metadata, transfer matrices, graph components, and safe partition exits.
%
% ilya.kuprov@weizmann.ac.il

function result=test_dynamic_metadata_partition_suite()

% State the metadata and partition target of the test
result=new_test_result('kernel/dynamic_metadata_partition_suite',...
                       'Metadata, hashing, and partition utilities',...
                       'small metadata and partition helpers must preserve stable identity and exact graph or algebraic behaviour.');

% Check parallel-state metadata on the MATLAB client
pool_count=poolsize();
result=test_true(result,'poolsize client scalar',...
                 isnumeric(pool_count)&&isscalar(pool_count)&&(pool_count>=0)&&...
                 (mod(pool_count,1)==0),...
                 'poolsize must return a non-negative integer worker count on the client');
result=test_true(result,'isworkernode client',~isworkernode(),...
                 'the validation driver runs on the MATLAB client rather than a parallel worker');

% Check MD5 hash stability and object-type sensitivity
hash_a=md5_hash({[1 2 3],'abc'});
hash_b=md5_hash({[1 2 3],'abc'});
hash_c=md5_hash({[1 2 4],'abc'});
result=test_true(result,'md5_hash stable hex',strcmp(hash_a,hash_b)&&numel(hash_a)==32&&all(isstrprop(hash_a,'xdigit')),...
                 'identical MATLAB objects must produce identical 32-character hexadecimal MD5 hashes');
result=test_true(result,'md5_hash data sensitivity',~strcmp(hash_a,hash_c),...
                 'changing the serialised object contents must change the MD5 hash');
result=test_true(result,'md5_hash storage sensitivity',~strcmp(md5_hash(eye(2)),md5_hash(speye(2))),...
                 'full and sparse matrices are distinct MATLAB objects and must hash differently');

% Check stable duplicate-row removal through hash-table identity
A=sparse([1 0 2;1 0 2;0 3 0;1 0 2;4 0 0]);
result=test_close(result,'unihash stable rows',unihash(A),sparse([1 0 2;0 3 0;4 0 0]),1e-15,1e-15,...
                  'unihash must keep the first occurrence of each unique sparse row in stable order');

% Check least-squares transfer matrix recovery from overdetermined samples
T_ref=[2 1;0 -1];
amp_inps=[1 0 1;0 1 1];
amp_outs=T_ref*amp_inps;
result=test_close(result,'transfermat exact recovery',transfermat(amp_inps,amp_outs),T_ref,1e-14,1e-14,...
                  'linearly complete input-output samples must recover the exact linear transfer matrix');

% Check strongly connected components on a two-component directed graph
G=logical([1 1 0;1 1 0;0 0 1]);
sci=scomponents(G);
result=test_true(result,'scomponents component partition',sci(1)==sci(2)&&sci(3)~=sci(1)&&numel(unique(sci))==2,...
                 'nodes one and two are mutually reachable, and node three is a separate component');

% Check path tracing disabled exit without graph partition work
spin_system.sys.output='hush';
spin_system.sys.disable={'pt'};
projectors=path_trace(spin_system,speye(3),[]);
result=test_true(result,'path_trace disabled projector',iscell(projectors)&&isscalar(projectors)&&isequal(projectors{1},1),...
                 'disabled path tracing must return a unit projector placeholder');

% Check zero-track elimination disabled exit without Krylov propagation
spin_system.bas.formalism='sphten-liouv';
spin_system.sys.disable={'zte'};
projector=zte(spin_system,speye(3),[1;0;0]);
result=test_true(result,'zte disabled projector',isequal(projector,1),...
                 'disabled zero-track elimination must return a unit projector placeholder');

end


