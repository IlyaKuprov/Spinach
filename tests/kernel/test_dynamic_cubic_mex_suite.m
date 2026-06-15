% Tests the cubic-polynomial MEX helper used by eigenfields(). Syntax:
%
%                    result=test_dynamic_cubic_mex_suite()
%
% Outputs:
%
%     result  - regression test result with explanatory messages
%
% The test checks cubic roots, degenerate lower-order polynomials,
% repeated roots, endpoint roots, extreme coefficient scaling, and
% derivative-root use cases against explicit references.
%
% ilya.kuprov@weizmann.ac.il

function result=test_dynamic_cubic_mex_suite()

% Announce the test target
fprintf('TESTING: Cubic polynomial root MEX helper\n');

% State the utility target of the test
result=new_test_result('kernel/dynamic_cubic_mex_suite',...
                       'Cubic polynomial root MEX helper',...
                       'eigenfields cubic root helper must match analytical and Matlab references.');

% Set the production root tolerance
root_tol=sqrt(eps);

% Check three roots, including endpoints
result=test_close(result,'cubic_roots three endpoint roots',...
                  cubic_roots([1 -1.5 0.5 0],root_tol),[0 0.5 1],...
                  1e-12,1e-12,...
                  'x*(x-1/2)*(x-1) must return all three unit-interval roots');

% Check a triple root
result=test_close(result,'cubic_roots triple root',...
                  cubic_roots([1 -1.5 0.75 -0.125],root_tol),0.5,...
                  1e-9,1e-12,...
                  'a cubic with a triple root should return one merged root');

% Check a double root plus an endpoint root
result=test_close(result,'cubic_roots double root',...
                  cubic_roots([1 -1.4 0.49 0],root_tol),[0 0.7],...
                  1e-9,1e-12,...
                  'x*(x-0.7)^2 should return the endpoint and merged double root');

% Check quadratic and linear degeneracies
result=test_close(result,'cubic_roots quadratic degeneracy',...
                  cubic_roots([0 1 -1 0],root_tol),[0 1],...
                  1e-12,1e-12,...
                  'leading zero cubic coefficient should reduce to the quadratic branch');
result=test_close(result,'cubic_roots linear degeneracy',...
                  cubic_roots([0 0 2 -1],root_tol),0.5,...
                  1e-12,1e-12,...
                  'leading zero cubic and quadratic coefficients should reduce to the linear branch');

% Check constant and zero polynomials
result=test_true(result,'cubic_roots constant polynomial',isempty(cubic_roots([0 0 0 1],root_tol)),...
                 'a non-zero constant polynomial has no roots');
result=test_true(result,'cubic_roots zero polynomial',isempty(cubic_roots([0 0 0 0],root_tol)),...
                 'the identically zero polynomial is ignored by the eigenfields root filter');

% Check extreme coefficient scaling
result=test_close(result,'cubic_roots large coefficient scale',...
                  cubic_roots([1e300 -1.5e300 0.5e300 0],root_tol),[0 0.5 1],...
                  1e-12,1e-12,...
                  'normalisation should prevent overflow on very large coefficients');
result=test_close(result,'cubic_roots small coefficient scale',...
                  cubic_roots([1e-300 -1.5e-300 0.5e-300 0],root_tol),[0 0.5 1],...
                  1e-12,1e-12,...
                  'normalisation should prevent underflow on very small coefficients');

% Check derivative-root use case
turn_ref=sort((3+[-1 1]*sqrt(3))/6);
result=test_close(result,'cubic_roots derivative quadratic',...
                  cubic_roots([0 3 -3 0.5],root_tol),turn_ref,...
                  1e-12,1e-12,...
                  'the derivative quadratic should return both turning points');

% Check a random set of well separated roots
rng_state=rng;
rng_cleanup=onCleanup(@()rng(rng_state));
rng(1);
rand_ok=true; rand_err=0;
for n=1:200
    test_roots=sort(0.05+0.90*rand(1,3));
    if min(diff(test_roots))<1e-3
        continue;
    end
    obs_roots=cubic_roots(poly(test_roots),root_tol);
    if numel(obs_roots)~=3
        rand_ok=false;
        break;
    end
    rand_err=max(rand_err,max(abs(obs_roots-test_roots)));
end



result=test_true(result,'cubic_roots random cubic set',rand_ok&&(rand_err<1e-10),...
                 'well-separated random cubic roots should round-trip through polynomial coefficients');
clear('rng_cleanup');

end
