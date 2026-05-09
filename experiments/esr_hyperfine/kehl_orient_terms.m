% Kehl tensor projections for a selected orientation. Syntax:
%
%      terms=kehl_orient_terms(parameters,euler_angles)
%
% Parameters:
%
%   parameters       - Kehl ENDOR context parameter structure.
%   euler_angles     - Kehl orientation angles stored by the context.
%
% Outputs:
%
%   terms            - structure containing hyperfine, NQI, CS, and dipolar projections.
%
% February 2024 A. Kehl (akehl@gwdg.de)
% May 2026 Spinach integration
%
% <https://spindynamics.org/wiki/index.php?title=kehl_orient_terms.m>

function terms=kehl_orient_terms(parameters,euler_angles)

    % Check consistency
    grumble(parameters,euler_angles);

    % Convert the stored Kehl angles into the legacy direction matrix
    theta=-euler_angles(2);
    phi=-euler_angles(3);
    R=zeros(3,3);
    R(1,1)=cos(theta)*cos(phi);
    R(1,2)=cos(theta)*sin(phi);
    R(1,3)=-sin(theta);
    R(2,1)=-sin(phi);
    R(2,2)=cos(phi);
    R(3,1)=sin(theta)*cos(phi);
    R(3,2)=sin(theta)*sin(phi);
    R(3,3)=cos(theta);

    % Initialise orientation-dependent ENDOR terms
    n_endor=parameters.n_endor;
    terms.hf_zz=zeros(1,n_endor);
    terms.hf_zy=zeros(1,n_endor);
    terms.hf_zx=zeros(1,n_endor);
    terms.nqi=zeros(n_endor,3,3);
    terms.nqi_zz=zeros(1,n_endor);
    terms.cs_zz=zeros(1,n_endor);

    % Project hyperfine, quadrupolar, and chemical-shift tensors
    for n=1:n_endor
        A=R*parameters.hfc_matrix(3*n-2:3*n,:)*R';
        terms.hf_zz(n)=A(3,3);
        if parameters.Bterm
            terms.hf_zy(n)=A(3,2);
            terms.hf_zx(n)=A(3,1);
        end
        if parameters.nqi_active
            Q=R*parameters.nqi_matrix(3*n-2:3*n,:)*R';
            terms.nqi_zz(n)=Q(3,3);
            if parameters.Bterm
                terms.nqi(n,:,:)=Q;
            end
        end
        if parameters.cs_active
            CS=R*parameters.cs_matrix(3*n-2:3*n,:)*R';
            terms.cs_zz(n)=CS(3,3);
        end
    end

    % Project nuclear dipolar tensors
    n_dip=size(parameters.dipolar_matrix,1)/3;
    terms.dip_zz=zeros(1,n_dip);
    if parameters.dipolar_active
        for n=1:n_dip
            D=R*parameters.dipolar_matrix(3*n-2:3*n,:)*R';
            terms.dip_zz(n)=D(3,3);
        end
    end

end

% Consistency enforcement
function grumble(parameters,euler_angles)
    if ~isstruct(parameters)
        error('parameters must be a structure.');
    end
    if (~isfield(parameters,'n_endor'))||(~isnumeric(parameters.n_endor))||...
            (~isscalar(parameters.n_endor))||(parameters.n_endor<1)
        error('parameters.n_endor must be a positive integer.');
    end
    if (~isfield(parameters,'hfc_matrix'))||(~isnumeric(parameters.hfc_matrix))||...
            (size(parameters.hfc_matrix,1)~=3*parameters.n_endor)||...
            (size(parameters.hfc_matrix,2)~=3)
        error('parameters.hfc_matrix has incompatible dimensions.');
    end
    if (~isfield(parameters,'nqi_matrix'))||(~isnumeric(parameters.nqi_matrix))||...
            (size(parameters.nqi_matrix,1)~=3*parameters.n_endor)||...
            (size(parameters.nqi_matrix,2)~=3)
        error('parameters.nqi_matrix has incompatible dimensions.');
    end
    if isfield(parameters,'cs_active')&&parameters.cs_active
        if (~isfield(parameters,'cs_matrix'))||(~isnumeric(parameters.cs_matrix))||...
                (size(parameters.cs_matrix,1)~=3*parameters.n_endor)||...
                (size(parameters.cs_matrix,2)~=3)
            error('parameters.cs_matrix has incompatible dimensions.');
        end
    end
    if (~isfield(parameters,'dipolar_matrix'))||(~isnumeric(parameters.dipolar_matrix))||...
            (size(parameters.dipolar_matrix,2)~=3)||...
            (mod(size(parameters.dipolar_matrix,1),3)~=0)
        error('parameters.dipolar_matrix has incompatible dimensions.');
    end
    if (~isnumeric(euler_angles))||(~isreal(euler_angles))||...
            (numel(euler_angles)~=3)
        error('euler_angles must be a three-element real vector.');
    end
end

