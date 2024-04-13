% Compares the output of Spinach Clebsch-Gordan function with the 
% arbitrary precision results returned by Mathematica.
%
% i.kuprov@soton.ac.uk

function cg_test_small()

% Load the test file
load('cg_test_table_small.mat','cg_test_table');

% Loop over the lines of the table
for n=1:size(cg_test_table,1)
    
    % Compute the CG coefficient
    tic;
    cg=clebsch_gordan(cg_test_table(n,1),cg_test_table(n,2),...
                      cg_test_table(n,3),cg_test_table(n,4),...
                      cg_test_table(n,5),cg_test_table(n,6));
    disp([num2str(toc) ' seconds']);
                  
    % Compare with the Mathematica result
    difference=abs(cg_test_table(n,7)-cg);
    if abs(difference)<2*eps
        disp(cg_test_table(n,:));
        disp([num2str(difference) ' - PASS']);
        disp('--------------------------------------------------------------');
    else
        disp(cg_test_table(n,:));
        error([num2str(difference) ' - FAIL']);
    end

end

end

