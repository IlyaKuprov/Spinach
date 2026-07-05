% Tests permutation group database metadata. Syntax:
%
%                    result=test_perm_group_suite()
%
% Outputs:
%
%     result  - regression test result with explanatory messages
%
% The test checks the real-valued-character Abelian permutation subgroups
% against explicit element tables, explicit character tables, character
% orthogonality, closure, commutativity, and involution order.
%
% ilya.kuprov@weizmann.ac.il

function result=test_perm_group_suite()

% Announce the test target
fprintf('TESTING: Permutation group metadata\n');

% State the utility target of the test
result=new_test_result('kernel/perm_group_suite',...
                       'Permutation group metadata',...
                       'real-valued-character Abelian permutation subgroup tables must be exact and internally consistent.');

% Check the S3 Abelian subgroup table
G=perm_group('S3A');
result=test_close(result,'perm_group S3A elements',G.elements,[1 2 3;2 1 3],1e-15,1e-15,...
                  'S3A is represented by a two-element transposition subgroup');
result=test_close(result,'perm_group S3A characters',G.class_characters,[1 1;1 -1],1e-15,1e-15,...
                  'S3A has the real character table of C2');

% Check the S5 Abelian subgroup table
G=perm_group('S5A');
result=test_close(result,'perm_group S5A elements',G.elements,[1 2 3 4 5;2 1 4 3 5;3 4 1 2 5;4 3 2 1 5],1e-15,1e-15,...
                  'S5A embeds the S4A Klein four subgroup with the fifth point fixed');
result=test_close(result,'perm_group S5A characters',G.class_characters,[1 1 1 1;1 1 -1 -1;1 -1 1 -1;1 -1 -1 1],1e-15,1e-15,...
                  'S5A has the real character table of C2 x C2');

% Check the S6 Abelian subgroup table
G=perm_group('S6A');
elements_ref=[1 2 3 4 5 6;2 1 4 3 5 6;3 4 1 2 5 6;4 3 2 1 5 6;...
              1 2 3 4 6 5;2 1 4 3 6 5;3 4 1 2 6 5;4 3 2 1 6 5];
chars_ref=[1 1 1 1 1 1 1 1;...
           1 1 -1 -1 1 1 -1 -1;...
           1 -1 1 -1 1 -1 1 -1;...
           1 -1 -1 1 1 -1 -1 1;...
           1 1 1 1 -1 -1 -1 -1;...
           1 1 -1 -1 -1 -1 1 1;...
           1 -1 1 -1 -1 1 -1 1;...
           1 -1 -1 1 -1 1 1 -1];
result=test_close(result,'perm_group S6A elements',G.elements,elements_ref,1e-15,1e-15,...
                  'S6A is represented by the direct product of S4A and S2');
result=test_close(result,'perm_group S6A characters',G.class_characters,chars_ref,1e-15,1e-15,...
                  'S6A has the real character table of C2 x C2 x C2');

% Check structural consistency of all real-valued-character Abelian options
group_names={'S3A','S4A','S5A','S6A'};
group_degrees=[3 4 5 6];
group_orders=[2 4 4 8];
for n=1:numel(group_names)

    % Get the group under test
    G=perm_group(group_names{n});

    % Check scalar group metadata
    result=test_true(result,[group_names{n} ' scalar metadata'],...
                     (G.order==group_orders(n))&&(G.nclasses==group_orders(n))&&(G.n_irreps==group_orders(n)),...
                     'Abelian real-character subgroup metadata must match the group order');

    % Check class and irrep dimensions
    result=test_close(result,[group_names{n} ' singleton classes'],G.class_sizes,ones(1,group_orders(n)),1e-15,1e-15,...
                      'Abelian subgroup conjugacy classes are singletons');
    result=test_close(result,[group_names{n} ' irrep dimensions'],G.irrep_dims,ones(1,group_orders(n)),1e-15,1e-15,...
                      'Abelian subgroup irreducible representations are one-dimensional');

    % Check character reality and orthogonality
    result=test_true(result,[group_names{n} ' real characters'],...
                     isreal(G.class_characters)&&all(abs(G.class_characters(:))==1),...
                     'real-valued-character Abelian subgroup characters must be signs');
    result=test_close(result,[group_names{n} ' character orthogonality'],...
                      G.class_characters*G.class_characters.',group_orders(n)*eye(group_orders(n)),1e-15,1e-15,...
                      'Abelian subgroup character rows must be orthogonal');

    % Initialise multiplication checks
    group_closed=true; group_comm=true; group_invol=true; group_perm=true;
    identity=1:group_degrees(n);

    % Check every permutation and every product
    for k=1:G.order
        element=G.elements(k,:);
        group_perm=group_perm&&isequal(sort(element),identity);
        group_invol=group_invol&&isequal(element(element),identity);
        for m=1:G.order
            prod_left=G.elements(k,G.elements(m,:));
            prod_right=G.elements(m,G.elements(k,:));
            group_closed=group_closed&&ismember(prod_left,G.elements,'rows');
            group_comm=group_comm&&isequal(prod_left,prod_right);
        end
    end

    % Check permutation-table structure
    result=test_true(result,[group_names{n} ' permutation rows'],group_perm,...
                     'all group elements must be valid permutation rows');
    result=test_true(result,[group_names{n} ' closure'],group_closed,...
                     'the listed permutation rows must be closed under composition');
    result=test_true(result,[group_names{n} ' commutativity'],group_comm,...
                     'the listed permutation rows must commute');
    result=test_true(result,[group_names{n} ' involutions'],group_invol,...
                     'all real-character Abelian subgroup elements must square to identity');

end


end
