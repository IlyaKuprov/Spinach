% Commutation and product relations for Pauli and
% central transition operators.
%
% ilya.kuprov@weizmann.ac.il

function commutation_6()

% Accuracy threshold
tol=1e-10;

% Test SU(2) relations for Pauli operators
for mult=[2 3 4 7]

    % Build Pauli operators
    S=pauli(mult);

    % Check angular momentum commutators
    err_xy=norm(comm(S.x,S.y)-1i*S.z,'fro');
    err_yz=norm(comm(S.y,S.z)-1i*S.x,'fro');
    err_zx=norm(comm(S.z,S.x)-1i*S.y,'fro');

    % Check ladder operator definitions
    err_p=norm(S.p-(S.x+1i*S.y),'fro');
    err_m=norm(S.m-(S.x-1i*S.y),'fro');

    % Report Pauli operator failures
    if max([err_xy err_yz err_zx err_p err_m])>tol
        error('Pauli operator commutation test FAILED.');
    end

end
disp('Pauli operator commutation test PASSED.');

% Test SU(2) relations for central transition operators
for mult=[4 6 8]

    % Build central transition operators
    CTx=centrans(mult,'x');
    CTy=centrans(mult,'y');
    CTz=centrans(mult,'z');
    CTp=centrans(mult,'+');
    CTm=centrans(mult,'-');

    % Build central transition support projector
    P=zeros(mult,mult);
    P(mult/2,mult/2)=1;
    P(mult/2+1,mult/2+1)=1;

    % Check commutators and product identities
    err_comm=norm(comm(CTx,CTy)-1i*CTz,'fro');
    err_lad_1=norm(comm(CTz,CTp)-CTp,'fro');
    err_lad_2=norm(comm(CTz,CTm)+CTm,'fro');
    err_lad_3=norm(comm(CTp,CTm)-2*CTz,'fro');
    err_prod=norm(CTp*CTm-(CTz+P/2),'fro');

    % Report central transition failures
    if max([err_comm err_lad_1 err_lad_2 err_lad_3 err_prod])>tol
        error('CT commutation test FAILED.');
    end

end
disp('CT commutation test PASSED.');

% Test central transition IST expansions
for mult=[4 6 8]
    for type={'x','y','z','+','-'}

        % Build central transition operator
        C=centrans(mult,type{1});

        % Obtain IST expansion.
        [states,coeffs]=ct2ist(mult,type{1});
        T=irr_sph_ten(mult);

        % Reconstruct the operator from IST basis terms
        C_rec=zeros(mult,mult,'like',1i);
        for n=1:numel(states)
            C_rec=C_rec+coeffs(n)*T{states(n)+1};
        end

        % Report IST expansion failures
        if norm(C-C_rec,'fro')>tol
            error('CT IST expansion test FAILED.');
        end

    end
end
disp('CT IST expansion test PASSED.');

end

