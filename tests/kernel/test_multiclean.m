% Tests fused panel multiplication against the existing clean_up. Syntax:
%
%                         result=test_multiclean()
%
% Outputs:
%
%     result  - regression result with per-case numerical,
%               precision, and shape checks
%
% Uses exactly representable dyadic inputs to isolate clean-up semantics
% from legitimate last-bit differences between blocked BLAS products.
%
% ilya.kuprov@weizmann.ac.il

function result=test_multiclean()

% State the test target
result=new_test_result('kernel/multiclean','Fused multiplication and clean-up',...
                       'Panel products agree with clean_up within numerical roundoff.');

% Exercise the supplied global storage policy
spin_system.sys.disable={}; spin_system.sys.parallel={'processes',1};
spin_system.tols.small_matrix=200;
spin_system.tols.dense_matrix=0.15;



% Exercise every sparse/full combination in both floating precisions
A=reshape(mod(1:35,7)-3,5,7)/8;
B=reshape(mod(1:77,9)-4,7,11)/8;
for precision={'double','single'}
    for complex_case=0:1
        for left_sparse=0:1
            for right_sparse=0:1
                left=cast(A,precision{1}); right=cast(B,precision{1});
                if complex_case, left=left+1i*left/2; right=right-1i*right; end
                if left_sparse, left=sparse(left); end
                if right_sparse, right=sparse(right); end
                for tol=[0 1/8 10]
                    label=sprintf('%s complex=%d sparse=%d%d tol=%g',...
                                  precision{1},complex_case,left_sparse,right_sparse,tol);
                    reference=reference_clean(spin_system,left*right,tol);
                    result=compare(result,label,multiclean(spin_system,left,right,tol),reference);
                end
            end
        end
    end
end

% Exercise density, dimension, and zero-matrix representation boundaries
for dimension=[199 200 201]
    for occupied=[0 1 30 31 201]
        column_count=min(occupied,dimension);
        A=speye(dimension);
        B=sparse([ones(dimension,column_count) zeros(dimension,201-column_count)]);
        A=sparse(A); B=sparse(B);
        for dense_input=0:1
            if dense_input, A=full(A); B=full(B); end
            label=sprintf('size=%d occupied=%d full=%d',dimension,occupied,dense_input);
            reference=reference_clean(spin_system,A*B,1/8);
            result=compare(result,label,multiclean(spin_system,A,B,1/8),reference);
        end
    end
end

% Exercise exactly the 0.15 density boundary on a sufficiently large matrix
A=speye(220); B=sparse([ones(220,33) zeros(220,187)]);
result=compare(result,'sparse exact density',multiclean(spin_system,A,B,1/8),...
               reference_clean(spin_system,A*B,1/8));
result=compare(result,'full exact density',multiclean(spin_system,full(A),full(B),1/8),...
               reference_clean(spin_system,full(A)*full(B),1/8));

% Clean-up rounds components, including ties, rather than dropping by magnitude
A=[-3/16 -1/16 1/16 3/16 1/32+1i/16];
result=compare(result,'rounding ties and complex components',multiclean(spin_system,1,A,1/8),...
               reference_clean(spin_system,A,1/8));

% Repeating clean_up is not generally idempotent in floating point
A=-302204600072882.88;
result=compare(result,'exactly one rounding pass',multiclean(spin_system,A,1,0.1),...
               reference_clean(spin_system,A,0.1));

% Cancellation is completed before any rounding is applied
A=[1 1]; B=[3/64 3/64;3/64 -3/64];
result=compare(result,'complete inner products',multiclean(spin_system,A,B,1/8),...
               reference_clean(spin_system,A*B,1/8));

% Compare ordinary non-dyadic products in every sparse/full combination
rng(7);
for sparse_a=[false true]
    for sparse_b=[false true]
        A=randn(23,17)+1i*randn(23,17); B=randn(17,31)+1i*randn(17,31);
        if sparse_a, A=sparse(A); end
        if sparse_b, B=sparse(B); end
        for tol=[0 1e-8]
            result=compare(result,'non-dyadic product',...
                           multiclean(spin_system,A,B,tol),...
                           reference_clean(spin_system,A*B,tol));
        end
    end
end

% Force the panel path by declaring an enormous worker count
rng(11); A=sprandn(200,120,0.05); B=sprandn(120,170,0.05);
many_workers=spin_system; many_workers.sys.parallel={'processes',1e12};
for tol=[0 1e-3]
    result=compare(result,'forced panels',multiclean(many_workers,A,B,tol),...
                   reference_clean(spin_system,A*B,tol));
    result=compare(result,'forced panels full',multiclean(many_workers,full(A),B,tol),...
                   reference_clean(spin_system,full(A)*B,tol));
end

% Exercise scalar multiplication, vectors, and zero inner dimensions
cases={{2,sparse(A)},{sparse(A),2},{single(2),single(B)},...
       {zeros(0,3),zeros(3,4)},{zeros(4,3),zeros(3,0)},...
       {zeros(4,0),zeros(0,5)},{sparse(4,0),sparse(0,5)},...
       {zeros(0,0),2},{2,zeros(0,0)},{(1:4)',1:7},...
       {single(eye(3)),double(eye(3))},...
       {sparse([1 0 1]),sparse(reshape(mod(1:900,7),3,300))},...
       {sparse(reshape(mod(1:900,7),300,3)),sparse([1;0;1])}};
for case_index=1:numel(cases)
    A=cases{case_index}{1}; B=cases{case_index}{2};
    result=compare(result,['shape case ' num2str(case_index)],...
                   multiclean(spin_system,A,B,1/8),reference_clean(spin_system,A*B,1/8));
end

% A single tolerance can change precision through existing clean_up arithmetic
for sparse_case=0:1
    A=eye(220); if sparse_case, A=sparse(A); end
    for tol=single([0 1/8])
        result=compare(result,'single tolerance precision',multiclean(spin_system,A,A,tol),...
                       reference_clean(spin_system,A*A,tol));
    end
end

% Preserve nonfinite arithmetic rather than filtering exceptional entries
A=[Inf NaN;1 -Inf]; B=eye(2);
result=compare(result,'nonfinite entries',multiclean(spin_system,A,B,1/8),...
               reference_clean(spin_system,A*B,1/8));

% Repeated cleanup of large non-dyadic entries may differ by roundoff
for sparse_case=0:1
    A=[-302204600072882.88 302204600072882.88];
    if sparse_case, A=sparse(A); end
    once=reference_clean(spin_system,A,0.1); twice=reference_clean(spin_system,once,0.1);
    result=test_true(result,'adversarial rounding is non-idempotent',...
                     ~isequaln(once,twice),'The case must distinguish repeated rounding.');
    result=compare(result,'single rounding of large entries',...
                   multiclean(spin_system,1,A,0.1),once);
end

% Disabled clean-up returns native BLAS arithmetic, not panel arithmetic
rng(17); A=randn(31,47); B=randn(47,53);
result=compare(result,'zero tolerance numerical agreement',multiclean(spin_system,A,B,0),A*B);

% Compare unrounded products within floating-point roundoff
rng(17); A=randn(31,47); B=randn(47,53);
spin_system.sys.disable={'clean-up'};
result=compare(result,'disabled cleanup numerical agreement',multiclean(spin_system,A,B,1),A*B);
spin_system.sys.disable={};

% Global density differs from each panel's local density
A=speye(240); B=sparse(240,40); B(:,1:4)=1;
for density=[0.05 0.15 0.9]
    for small_matrix=[0 200]
        spin_system.tols.dense_matrix=density;
        spin_system.tols.small_matrix=small_matrix;
        A=sparse(A); B=sparse(B);
        for dense_input=0:1
            if dense_input, A=full(A); B=full(B); end
            label=sprintf('global density=%g size=%d full=%d',density,small_matrix,dense_input);
            result=compare(result,label,multiclean(spin_system,A,B,0.125),...
                           clean_up(spin_system,A*B,0.125));
        end
    end
end
spin_system.tols.dense_matrix=0.15;
spin_system.tols.small_matrix=200;

% Custom decoupling retains its unfused class-aware cleanup
sys.isotopes={'1H','13C'}; sys.output='hush'; sys.magnet=9.4;
sys.disable={'hygiene'}; sys.parallel={'processes',1};
inter.zeeman.scalar={0 0};
for formalism={'sphten-liouv','zeeman-liouv','zeeman-hilb'}
    bas.formalism=formalism{1}; bas.approximation={'none'};
    custom_system=basis(create(sys,inter),bas);
    A=operator(custom_system,'Lz','1H')+operator(custom_system,'Lz','13C');
    if ~strcmp(formalism{1},'zeeman-hilb')
        for disabled=[false true]
            custom_system.sys.disable={'hygiene'};
            if disabled, custom_system.sys.disable{end+1}='clean-up'; end
            numeric=decouple(custom_system,A,[],{'1H'});
            custom=decouple(custom_system,polyadic({{A}}),[],{'1H'});
            label=sprintf('custom fallback %s disabled=%d',formalism{1},disabled);
            result=compare(result,label,full(custom),full(numeric));
        end
    end
end

% Preserve custom multiplication and cleanup for either operand
custom=polyadic({{sparse(diag([1 2 3 4])/8)}});
operands={{custom,eye(4)},{eye(4),custom},{custom,custom}};
for pair=1:numel(operands)
    for disabled=[false true]
        custom_system=spin_system; custom_system.sys.disable={};
        if disabled, custom_system.sys.disable={'clean-up'}; end
        for tol=[0 1/8]
            left=operands{pair}{1}; right=operands{pair}{2};
            expected=clean_up(custom_system,left*right,tol);
            actual=multiclean(custom_system,left,right,tol);
            result=test_true(result,'custom operand values',...
                             isequaln(full(actual),full(expected)),...
                             'Custom objects must retain native multiplication and cleanup.');
            result=test_true(result,'custom operand representation',...
                             strcmp(class(actual),class(expected)),...
                             'Custom representation must not be materialised by multiclean.');
        end
    end
end

% Reject malformed inputs without accepting ambiguous tolerance policies
bad_cases={{ones(2),ones(2),NaN},{ones(2),ones(2),Inf},...
           {ones(2),ones(2),int32(0)},{ones(2),ones(2),uint64(1)},...
           {ones(2),ones(3),1},{ones(2),ones(2),-1},...
           {ones(2),ones(2),1i},{ones(2),ones(2),[1 2]},...
           {ones(2,2,2),1,1}};
for case_index=1:numel(bad_cases)
    rejected=false;
    try
        multiclean(spin_system,bad_cases{case_index}{:});
    catch
        rejected=true;
    end
    result=test_true(result,['invalid case ' num2str(case_index)],rejected,...
                     'Invalid arguments must be rejected.');
end

end

% Compare numerical values to native multiplication followed by cleanup
function result=compare(result,label,observed,reference)

% Check consistency
grumble(result,label,observed,reference);

% Allow accumulated floating-point roundoff relative to the reference
result=test_true(result,[label ' values'],isequaln(full(observed),full(reference))||...
                 norm(full(observed-reference),'fro')<=...
                 max(size(reference))*eps(class(reference))*max(1,norm(full(reference),'fro')),...
                 'Panel clean-up must agree with clean_up of the complete product.');
result=test_true(result,[label ' precision'],strcmp(class(observed),class(reference)),...
                 'The multiplication precision must be preserved.');

end

% Consistency enforcement
function grumble(result,label,observed,reference)
if (~isstruct(result))||(~ischar(label))||(~isrow(label))
    error('a test result and a character row label are required.');
end
if (~isnumeric(observed))||(~isnumeric(reference))
    error('observed and reference must be numeric arrays.');
end
end


% Reference rounding and global storage policy
function C=reference_clean(spin_system,C,tol)

C=clean_up(spin_system,C,tol);

end


