% Exhaustive product-state sampling and clock-transition regression tests.
% Syntax: result=test_cce_sampling()
%
% ilya.kuprov@weizmann.ac.il

function result=test_cce_sampling()

% Build a central avoided crossing with two longitudinal bath fields
result=new_test_result('kernel/cce_sampling','CCE product-state sampling',...
                       'Compare weighted samples, exact mixtures and clock-transition exterior fields.');
[spin_system,p]=fixture(); model=cce_prepare(spin_system,p);
p.bath_state=struct('mode','sampled','samples',1,'seed',17,'basis','uniform',...
                    'probabilities',{{[0.7;0.3],[0.2;0.8]}});
ensemble=cce_states(model,p); schedule=cce_sequence(model,p);
weights=[0.14 0.56 0.06 0.24]; sampled=zeros(2,2,numel(p.times));
for a=1:2
    for b=1:2
        k=2*(a-1)+b; context=ensemble.contexts(1); assignment=[a b];
        for j=1:2
            state=zeros(2); state(assignment(j),assignment(j))=1;
            context.bath_rho{j}=state; context.expectations(j,:)=[0 0 1.5-assignment(j)];
        end
        raw=cce_cluster(model,struct('id',1,'sites',[2 3]),context,schedule,p);
        sampled=sampled+weights(k)*raw.density;
    end
end
p.bath_state=struct('mode','explicit_product','rho',{{diag([0.7 0.3]),diag([0.2 0.8])}});
ensemble=cce_states(model,p); exact=cce_cluster(model,struct('id',1,'sites',[2 3]),ensemble.contexts(1),schedule,p);
mixture_error=norm(sampled(:)-exact.density(:));
result=test_close(result,'Exhaustive weighted samples equal exact product mixture',sampled,exact.density,1e-12,0,'Exhaustive weighted samples equal exact product mixture');

% Show that reducing samples and then averaging is not the reversed estimator
sample_products=[2 1/2].*[3 1/3]; average_products=mean(sample_products);
product_averages=mean([2 1/2])*mean([3 1/3]);
result=test_true(result,'Average of products differs from product of averages',abs(average_products-product_averages)>0.9,'Average of products differs from product of averages');

% Compare static exterior fields to exact quantum dynamics at an avoided crossing
context=ensemble.contexts(1); context.bath_rho{2}=diag([1 0]); context.expectations(2,:)=[0 0 0.5];
full=cce_cluster(model,struct('id',1,'sites',[2 3]),context,schedule,p);
p.meanfield='static'; sampled_field=cce_cluster(model,struct('id',1,'sites',2),context,schedule,p);
p.meanfield='none'; absent_field=cce_cluster(model,struct('id',1,'sites',2),context,schedule,p);
field_error=norm(sampled_field.density(:)-full.density(:));
field_effect=norm(absent_field.density(:)-full.density(:));
result=test_close(result,'Clock-transition static eigenstate exterior matches exact retained bath',sampled_field.density,full.density,1e-12,0,'Clock-transition static eigenstate exterior matches exact retained bath');
result=test_true(result,'Omitting sampled exterior field changes clock-transition dynamics',field_effect>0.01,'Omitting sampled exterior field changes clock-transition dynamics');

% Check the avoided-crossing energy is stationary to first order in longitudinal field
spin=pauli(2); x=spin.x; z=spin.z; delta=1e-4; gap=@(b) diff(sort(eig(5*x+b*z)));
derivative=(gap(delta)-gap(-delta))/(2*delta);
result=test_close(result,'Central gap has zero first-order longitudinal sensitivity',derivative,0,1e-12,0,'Central gap has zero first-order longitudinal sensitivity');
fprintf('sampling mixture error %.3g; estimator difference %.6g; clock field error %.3g; omitted-field effect %.6g\n',...
        mixture_error,average_products-product_averages,field_error,field_effect);

end

% Prepare the central transverse splitting using canonical lab-frame metadata
function [spin_system,p]=fixture()
sys.isotopes={'1H','1H','1H'}; sys.magnet=0; sys.disable={'hygiene'}; sys.parallel={'processes',1}; sys.output='hush';
inter.coupling.matrix=cell(3); inter.coupling.matrix{1,2}=diag([0 0 3])/(2*pi);
inter.coupling.matrix{1,3}=diag([0 0 4])/(2*pi); spin_system=create(sys,inter);
spin_system.inter.zeeman.matrix{1}=[0 0 5;0 0 0;5 0 0];
p.central_spins=1; p.bath_spins=[2 3]; p.method='gcce'; p.expansion='standard';
p.quantity='density'; p.times=[0 0.2 0.4]; p.order=2; p.central.rho=diag([1 0]);
p.assumptions='labframe'; p.orientation=[0 0 0]; p.interaction_cutoff=Inf;
p.graph=struct('mode','complete','score','coupling_norm'); p.meanfield='none'; p.sequence=[];
p.numerics=struct('atol',1e-9,'rtol',1e-9,'max_memory_bytes',2^26);
end

