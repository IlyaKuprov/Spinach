% Declared-window envelope fitting and invalid-data diagnostics.
%
% ilya.kuprov@weizmann.ac.il

function result=test_cce_analysis()

% Fit a complex coherence through its explicitly selected magnitude envelope
result=new_test_result('kernel/cce_analysis','CCE decay diagnostics',...
                       'Synthetic windows, fit assumptions and invalid/oscillatory refusal.');
t=linspace(0,2,101)'; envelope=0.8*exp(-(t/0.7).^1.6);
r=struct('times',t,'value',envelope.*exp(4i*t),'valid',true(size(t)));
p.analysis=struct('fit',true,'observable',1,'envelope','magnitude',...
                  'model','stretched_exponential','interval',[0.1 1.8],...
                  'weights',ones(size(t)),'confidence','linearized','confidence_level',0.95);
a=cce_analyse(r,p);
result=test_true(result,'C7 contract check',strcmp(a.status,'fitted'),'Declared envelope fits successfully.');
result=test_close(result,'C7 numerical check',[a.fit.amplitude a.fit.time_constant a.fit.beta],[0.8 0.7 1.6],...
                  1e-7,1e-7,'Synthetic stretched exponential parameters.');
result=test_true(result,'C7 contract check',isequal(a.value,r.value)&&isequal(a.valid,r.valid),'Raw curve and masks retained.');
result=test_true(result,'C7 contract check',isequal(size(a.fit.confidence_interval),[3 2]),'Declared confidence procedure recorded.');

% Confirm fit-off default and refusal to hide an invalid sample
idle=cce_analyse(r,struct()); result=test_true(result,'C7 contract check',strcmp(idle.status,'not_requested'),'Fitting is opt-in.');
r.valid(20)=false; invalid=cce_analyse(r,p);
result=test_true(result,'C7 contract check',strcmp(invalid.status,'invalid_window')&&~isempty(invalid.warnings),...
                 'Invalid samples are not silently removed from the fitting window.');

% Refuse an oscillatory real curve rather than inventing a relaxation time
r.valid(:)=true; r.value=exp(-t).*(1+0.3*cos(20*t)); p.analysis.envelope='real';
oscillatory=cce_analyse(r,p);
result=test_true(result,'C7 contract check',strcmp(oscillatory.status,'not_decay')&&~isempty(oscillatory.warnings),...
                 'Oscillatory windows require an explicit valid envelope.');
fprintf('Fit A %.10g T %.10g beta %.10g; weighted RMS %.3g\n',...
        a.fit.amplitude,a.fit.time_constant,a.fit.beta,a.fit.weighted_rms);

end


