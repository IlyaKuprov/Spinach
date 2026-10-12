% Tests parallel pool ValueStore snapshot saving and merging. Syntax:
%
%                         result=test_vstore()
%
% Outputs:
%
%     result - regression checks for saved values, retained keys,
%              overwritten keys, empty snapshots, and invalid files
%
% talos@spindynamics.org

function result=test_vstore()

% State the snapshot persistence contract
result=new_test_result('kernel/vstore','ValueStore snapshot merging',...
                      'loading a snapshot must add its values without clearing unrelated keys.');

% Create a real shared store through Spinach
sys.magnet=1; sys.isotopes={'1H'}; sys.output='hush';
sys.parallel={'processes',1}; sys.parprops={}; sys.disable={'hygiene'};
inter.zeeman.scalar={0}; create(sys,inter);
pool=gcp('nocreate'); store=pool.ValueStore;
file_name=[tempname '.mat'];
cleanup_obj=onCleanup(@()delete(file_name));

% Save heterogeneous values and check that saving leaves them intact
saved_keys=["matrix";"sparse";"nested";"shared"];
saved_vals={complex([1 2;3 4],[4 3;2 1]);speye(3);...
            struct('items',{{"text",[],uint64(17)}});"saved"};
put(store,saved_keys,saved_vals); save_vstore(file_name);
result=test_true(result,'save retains values',isequal(get(store,saved_keys),saved_vals),...
                 'saving must not modify the current store');

% Merge into a populated store with a conflicting key
remove(store,keys(store)); put(store,["current";"shared"],{42;"old"});
load_vstore(file_name);
result=test_true(result,'snapshot values restored',isequal(get(store,saved_keys),saved_vals),...
                 'dense, sparse, nested, and string values must round-trip exactly');
result=test_true(result,'sparse representation retained',issparse(store("sparse")),...
                 'MAT-file persistence must retain sparse storage');
result=test_true(result,'current-only key retained',isequal(store("current"),42),...
                 'loading must not remove a key absent from the snapshot');
result=test_true(result,'matching key overwritten',isequal(store("shared"),"saved"),...
                 'the saved value must replace the current value for the same key');
result=test_true(result,'merged key set',isequal(sort(keys(store)),sort([saved_keys;"current"])),...
                 'the merged store must contain the union of current and saved keys');

% Repeat the merge without changing its result
load_vstore(file_name);
result=test_true(result,'repeated merge',isequal(get(store,saved_keys),saved_vals)&&...
                 isequal(store("current"),42),...
                 'loading the same snapshot again must retain the merged values');

% Load an empty snapshot into a populated store
remove(store,keys(store)); save_vstore(file_name); put(store,"current",{42});
load_vstore(file_name);
result=test_true(result,'empty snapshot preserves store',isequal(keys(store),"current")&&...
                 isequal(store("current"),42),...
                 'an empty snapshot must leave the current store unchanged');

% Reject a malformed snapshot before altering the current store
key_set="invalid"; val_set=17; save(file_name,'key_set','val_set');
rejected=false;
try
    load_vstore(file_name);
catch err
    rejected=contains(err.message,'ValueStore snapshot from save_vstore');
end
result=test_true(result,'invalid snapshot rejected',rejected,...
                 'non-cell snapshot values must be rejected');
result=test_true(result,'invalid snapshot preserves store',isequal(keys(store),"current")&&...
                 isequal(store("current"),42),...
                 'snapshot format errors must not alter the store');

end


