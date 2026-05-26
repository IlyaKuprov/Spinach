% Test of matrix exponential differentiation routines. Analytical
% derivatives are compared to central finite differences.
%
% ilya.kuprov@weizmann.ac.il

function dirdiff_1()

% Formalisms to test
formalisms={'sphten-liouv','zeeman-liouv','zeeman-hilb'};

% Loop over formalisms
for n=1:numel(formalisms)

    % Get the Spinach object
    spin_system=dirdiff_test_system(formalisms{n});

    % Random Hamiltonian
    H=randn(5)+1i*randn(5); H=(H+H')/2;

    % Random direction operators
    A=randn(5)+1i*randn(5); A=(A+A')/20;
    B=randn(5)+1i*randn(5); B=(B+B')/20;

    % First derivative, numerical
    D_num=(propagator(spin_system,H+1e-3*A,1)-...
           propagator(spin_system,H-1e-3*A,1))/2e-3;

    % First derivative, analytical
    D_anl=dirdiff(spin_system,H,A,1,2);

    % Test the first derivative
    if norm(D_num-D_anl{2},2)/norm(D_num,2)>1e-5
        error([formalisms{n} ' first derivative test failed']);
    else
        disp([formalisms{n} ' first derivative test passed']);
    end

    % Second derivative, numerical
    D_num=(propagator(spin_system,H+1e-3*A+1e-3*B,1)-...
           propagator(spin_system,H+1e-3*A-1e-3*B,1)-...
           propagator(spin_system,H-1e-3*A+1e-3*B,1)+...
           propagator(spin_system,H-1e-3*A-1e-3*B,1))/4e-6;

    % Second derivative, analytical
    P=dirdiff(spin_system,H,{A,B},1,3);
    Q=dirdiff(spin_system,H,{B,A},1,3);
    D_anl=(P{3}+Q{3})/2;

    % Test the second derivative
    if norm(D_num-D_anl,2)/norm(D_num,2)>1e-3
        error([formalisms{n} ' second derivative test failed']);
    else
        disp([formalisms{n} ' second derivative test passed']);
    end

end

end
