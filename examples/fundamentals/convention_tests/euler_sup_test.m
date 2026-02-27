% Euler angle superposition tests.
%
% ilya.kuprov@weizmann.ac.il

function euler_sup_test()

% Identity composition test
ang=euler_sup([0 0 0],[0 0 0]);
if norm(euler2dcm(ang)-eye(3),1)>1e-12
    error('identity composition test failed.');
end

% Random stress test
for n=1:2000

    % Draw random Euler angles
    ang_one=8*pi*(rand(1,3)-0.5);
    ang_two=8*pi*(rand(1,3)-0.5);

    % Compose through Euler superposition utility
    ang_cmp=euler_sup(ang_one,ang_two);

    % Compose through direct matrix multiplication
    dcm_ref=euler2dcm(ang_two)*euler2dcm(ang_one);

    % Compare the two composite matrices
    if norm(euler2dcm(ang_cmp)-dcm_ref,1)>1e-3
        error('random stress test failed.');
    end

end

% Singular branch stress test A
for n=1:500

    % Draw random near-singular rotations
    ang_one=[8*pi*(rand-0.5) 1e-12*randn() 8*pi*(rand-0.5)];
    ang_two=[8*pi*(rand-0.5) 1e-12*randn() 8*pi*(rand-0.5)];

    % Compose through Euler superposition utility
    ang_cmp=euler_sup(ang_one,ang_two);

    % Compose through direct matrix multiplication
    dcm_ref=euler2dcm(ang_two)*euler2dcm(ang_one);

    % Compare the two composite matrices
    if norm(euler2dcm(ang_cmp)-dcm_ref,1)>1e-3
        error('singular beta~0 stress test failed.');
    end

end

% Singular branch stress test B
for n=1:500

    % Draw random near-singular rotations
    ang_one=[8*pi*(rand-0.5) pi+1e-12*randn() 8*pi*(rand-0.5)];
    ang_two=[8*pi*(rand-0.5) pi+1e-12*randn() 8*pi*(rand-0.5)];

    % Compose through Euler superposition utility
    ang_cmp=euler_sup(ang_one,ang_two);

    % Compose through direct matrix multiplication
    dcm_ref=euler2dcm(ang_two)*euler2dcm(ang_one);

    % Compare the two composite matrices
    if norm(euler2dcm(ang_cmp)-dcm_ref,1)>1e-3
        error('singular beta~pi stress test failed.');
    end

end

% Singular branch stress test C
for n=1:500

    % Draw same-phase rotations adding to beta=pi
    alp=2*pi*(rand-0.5);
    gam=2*pi*(rand-0.5);
    bet=pi*rand;
    ang_one=[alp bet gam];
    ang_two=[alp pi-bet gam];

    % Compose through Euler superposition utility
    ang_cmp=euler_sup(ang_one,ang_two);

    % Compose through direct matrix multiplication
    dcm_ref=euler2dcm(ang_two)*euler2dcm(ang_one);

    % Compare the two composite matrices
    if norm(euler2dcm(ang_cmp)-dcm_ref,1)>1e-3
        error('mrsimulator singular branch test failed.');
    end

end

% Successful completion message
disp('Euler angle superposition test PASSED.');

end

