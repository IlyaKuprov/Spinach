% Smoke and independent-reference checks for the CCE application examples.
% Syntax: result=test_cce_examples()
%
% ilya.kuprov@weizmann.ac.il

function result=test_cce_examples()

% Initialise the standard test record and preserve figure visibility
result=new_test_result('kernel/cce_examples','CCE application examples',...
                       'Runnable categories and independently evaluated small limits.');
visibility=get(groot,'defaultFigureVisible');
cleanup=onCleanup(@() set(groot,'defaultFigureVisible',visibility));
set(groot,'defaultFigureVisible','off');

% Enumerate examples explicitly so new unrelated files cannot enter the test
names={'cce_single_nucleus','cce_nv_diamond','cce_sip_donor',...
       'cce_sic_divacancy','gcce_sic_clock','cce_p1_electron_bath',...
       'pcce_electron_layer','gcce_votpp','cce_erbium_cawo4',...
       'gcce_highspin_siv','cce_hbn_quadrupoles','cce_polarized_bath',...
       'gcce_two_centres','gcce_finite_drive','cce_deer_bandwidth',...
       'cce_noise_spectrum','mecce_dissipative_bath',...
       'cce_methyl_molecular','cce_molecular_powder'};
for n=1:numel(names)

    % Execute each public example and retain the standard harness fields only
    started=tic; actual=feval(names{n});
    result=test_true(result,[names{n} ' valid'],all(actual.valid,'all'),...
                     'Every displayed time point has a valid expansion.');
    result=test_true(result,[names{n} ' finite'],all(isfinite(actual.value),'all'),...
                     'Every returned complex value is finite.');
    if isfield(actual,'reference_error')

        % Independent small-system targets use an absolute trace-error criterion
        result=test_close(result,[names{n} ' reference'],actual.reference_error,...
                          0,1e-7,0,'Analytic, dense expm, or pinned Python reference.');
    end
    if ~isempty(actual.density)

        % Full-order finite baths must preserve observed density physicality
        result=test_true(result,[names{n} ' trace'],...
                         max(actual.diagnostics.trace_error)<1e-8,...
                         'Reduced density matrices retain unit trace.');
        result=test_true(result,[names{n} ' positive'],...
                         min(actual.diagnostics.minimum_eigenvalue)>-1e-8,...
                         'No material negative population at full bath order.');
    end
    if strcmp(names{n},'cce_deer_bandwidth')

        % Excluding the pump resonance must remove its selected pulse
        parameters=actual.input_parameters;
        parameters.deer.detuning_center=2*pi*1e7;
        sys.magnet=0.01; sys.isotopes={'E','E','E'}; sys.output='hush';
        inter.zeeman.scalar={0,0,0};
        model=cce_prepare(create(sys,inter),parameters);
        schedule=cce_controls(model,[],[],[],parameters);
        result=test_true(result,'DEER example excludes pump resonance',...
                         isempty(schedule.selected_sites{2}),...
                         'The pump pulse must follow the documented frequency window.');
    end
    if strcmp(names{n},'cce_molecular_powder')

        % All public response arrays must describe the powder average
        result=test_close(result,'Powder response alias',actual.response,...
                          actual.value,0,0,...
                          'The response alias must equal the plotted powder signal.');
        result=test_close(result,'Powder final order response',...
                          actual.by_order(:,:,end),actual.value,1e-14,...
                          1e-14,'The final order must be orientation-averaged.');
        expected=zeros(size(actual.abs_error_bound));
        for k=1:numel(actual.orientation_weights)
            expected=expected+actual.orientation_weights(k)*...
                              actual.orientation_abs_error_bound(:,k);
        end
        result=test_close(result,'Powder averaged numerical bounds',...
                          actual.abs_error_bound,expected,0,0,...
                          'Positive quadrature weights also average absolute bounds.');
        result=test_close(result,'Powder diagnostic numerical bounds',...
                          actual.diagnostics.numerical_error,actual.abs_error_bound,...
                          0,0,...
                          'Diagnostics must be recomputed for the powder result.');
    end
    fprintf('Example %s completed in %.3f s.\n',names{n},toc(started));
    close all;
end

end


