% Tests waveform basis orthonormality. Syntax:
%
%                    result=test_wave_basis_orthonormality()
%
% Outputs:
%
%     result  - regression test result with explanatory messages
%
% The test checks that sine, cosine, and Legendre waveform bases returned
% by Spinach have orthonormal columns, as required by pulse optimisation.
%
% ilya.kuprov@weizmann.ac.il

function result=test_wave_basis_orthonormality()

% Announce the test target
fprintf('TESTING: Waveform basis orthonormality\n');

% State the numerical target of the test
result=new_test_result('kernel/wave_basis_orthonormality',...
                       'Waveform basis orthonormality',...
                       'pulse waveform basis columns must be orthonormal.');

% Check all supported basis families
basis_types={'sine_waves','cosine_waves','legendre'};
for n=1:numel(basis_types)

    % Build a small waveform basis
    B=wave_basis(basis_types{n},5,32);

    % Check orthonormality of columns
    result=test_close(result,[basis_types{n} ' Gram matrix'],B'*B,eye(5),1e-12,1e-12,...
                      'orthonormal columns give independent waveform coefficients');
end

end

