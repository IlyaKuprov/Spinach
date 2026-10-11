% Tests regression-runner output capture and graphics isolation. Syntax:
%
%                         result=test_runner_output()
%
% Outputs:
%
%     result - checks of passing, recorded-failure, and thrown-error paths

function result=test_runner_output()

% Protect caller graphics while constructing an existing user figure
cleaner=test_graphics(); %#ok<NASGU>
user_figure=figure('Visible','off');
hidden_figure=figure('Visible','off','HandleVisibility','off'); %#ok<NASGU>
user_axes=axes(user_figure);
line_handle=plot(user_axes,[1 2],[3 4]);
set(groot,'CurrentFigure',user_figure);
set(groot,'DefaultFigureVisible','on');
defaults=get(groot,'default');
figures=findall(groot,'Type','figure');

% Exercise the same execution path for success and both failure modes
result=new_test_result('runner/output_cleanup','Runner output and graphics cleanup',...
                       'test chatter and figures must not leak, and failures must survive capture.');
for mode=0:2
    entry=struct('id','fixture','name','Runner fixture',...
                 'function',@()runner_fixture(mode));
    [observed,transcript]=test_execute(entry);
    result=test_true(result,'captured chatter',contains(transcript,'fixture chatter'),...
                     'console output remains available to the failure reporter');
    result=test_true(result,'figure cleanup',isequal(findall(groot,'Type','figure'),figures),...
                     'ordinary and hidden test figures are deleted on every outcome');
    result=test_true(result,'graphics restoration',isequal(get(groot,'default'),defaults)&&...
                     isequal(get(groot,'CurrentFigure'),user_figure)&&isgraphics(line_handle),...
                     'root defaults, the current figure, and existing contents survive');
    result=test_true(result,'outcome',strcmp(observed.status,'PASS')==(mode==0),...
                     'recorded failures and thrown errors both fail the test');
    result=test_true(result,'partial checks',any(contains(observed.messages,'invisible')),...
                     'checks before a thrown error remain in the result');
    if mode==2
        result=test_true(result,'exception text',contains(observed.error,'fixture exception'),...
                         'the original thrown error is retained');
    end
end

% Verify cleanup also runs when the guarded function itself unwinds
try
    guarded_failure();
catch err
    result=test_true(result,'guard exception',strcmp(err.identifier,'Spinach:GuardFixture'),...
                     'the cleanup guard does not replace the original exception');
end
result=test_true(result,'unwinding cleanup',isequal(findall(groot,'Type','figure'),figures)&&...
                 isequal(get(groot,'default'),defaults)&&...
                 isequal(get(groot,'CurrentFigure'),user_figure),...
                 'an exception restores the graphics state without an explicit clear');

% Check the public reporter without recursively selecting this test
text=evalc('observed=run_tests(''pattern'',''kernel/pauli_spin_half_algebra'');');
result=test_true(result,'uniform progress',contains(text,'[1/1]')&&...
                 contains(text,'PASS')&&contains(text,'1 passed, 0 failed')&&...
                 ~contains(text,'TESTING:')&&~contains(text,'PASS:'),...
                 'normal reporting contains progress and summary but no test chatter');
text=evalc('observed=run_test(''kernel/pauli_spin_half_algebra'');');
result=test_true(result,'verbose checks',contains(text,'PASS:')&&...
                 ~contains(text,'TESTING:')&&strcmp(observed.status,'PASS'),...
                 'run_test adds recorded checks without exposing stray console output');

end

% Produce noise and graphics before each possible test outcome
function result=runner_fixture(mode)

% Record an invisible default figure and a hidden-handle figure
result=new_test_result('fixture','Runner fixture','Capture and cleanup probe.');
fprintf('fixture chatter\n');
fig=figure();
result=test_true(result,'invisible',strcmp(get(fig,'Visible'),'off'),...
                 'test figures inherit invisible visibility');
figure('HandleVisibility','off');
set(groot,'DefaultFigureColor',[0.1 0.2 0.3]);
if mode==1
    result=test_true(result,'fixture failure',false,'recorded failure survives');
elseif mode==2
    error('Spinach:RunnerFixture','fixture exception');
end

end

% Throw with a live graphics guard to exercise stack-unwinding cleanup
function guarded_failure()

% Leave both a figure and a changed default for the guard to restore
cleaner=test_graphics(); %#ok<NASGU>
figure('HandleVisibility','off');
set(groot,'DefaultFigureVisible','off');
error('Spinach:GuardFixture','guarded exception');

end


