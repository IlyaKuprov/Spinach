% Field-domain EPR orientation selection for Kehl ENDOR. Syntax:
%
%      EPR=kehl_ori_field(spin_system,parameters)
%
% Parameters:
%
%   spin_system      - Spinach spin system structure.
%   parameters       - Kehl ENDOR context parameter structure.
%
% Outputs:
%
%   EPR              - map containing selected field-domain orientations.
%
% February 2024 A. Kehl (akehl@gwdg.de)
% May 2026 Spinach integration
%
% <https://spindynamics.org/wiki/index.php?title=kehl_ori_field.m>

function EPR=kehl_ori_field(spin_system,parameters)

    % Check consistency
    grumble(spin_system,parameters);

    % Unpack context data
    constants=parameters.constants;
    paramsEPR=parameters.paramsEPR;
    paramsENDOR=parameters.paramsENDOR;
    Ntheta=parameters.Nang;
    Nphimax=parameters.Nang;
    g=parameters.g_matrix;
    A=parameters.hfc_matrix;
    Q=parameters.nqi_matrix;
    if parameters.epr_nuclei_active==true
        Ni_EPR=parameters.n_epr;
        A_EPR=parameters.epr_hfc_matrix;
    else
        Ni_EPR=0;
    end
    if parameters.cs_active==true
        CS=parameters.cs_matrix;
    end
    if parameters.dipolar_active==true
        D=parameters.dipolar_matrix;
    end
    n_endor=parameters.n_endor;
    g2=g*g;

    % initialize arrays
    tmp_epr=zeros(paramsEPR("Npts"),1);
    epr_amp=zeros(paramsEPR("Npts"),1);

    or=0;
    geff_sel=[];
    B_sel=[];
    HF_zz_sel=[];
    HF_zy_sel=[];
    HF_zx_sel=[];
    NQI_zz_sel=[];
    NQI_sel=[];
    CS_zz_sel=[];
    D_zz_sel=[];
    euler_sel=[];
    S_sel=[];
    offsets_sel=[];

    % for powder pattern
    if parameters.powder==true
        if isfield(parameters,'excite_width')
            W1=parameters.excite_width;
        else
            W1=parameters.pulse_width*constants("CONST1")*1e10;
        end
        nor=0;
        % Loop over orientations
        for ii=1:Ntheta
            theta=ii*pi/Ntheta;
            Nphi=round(sin(theta)*Nphimax)*1;
            for jj=1:Nphi
                phi=(jj-1)*pi*2/(Nphi);

                % Euler angles for Spinach orientation assembly
                euler_angles=[0 -theta -phi];

                % Direction cosine vector
                dc=[cos(phi)*sin(theta) sin(phi)*sin(theta) cos(theta)];

                % Effective g-value for a given theta, and phi combination
                geff=(dc*g2*dc')^.5;

                % effective B field for given theta and phi
                B=(parameters.field_t*parameters.g_iso)/geff(1);
                Beff=parameters.mw_freq_hz*constants('H')/(constants('MU_B')*geff);

                % HF ENDOR
                HF_zz=zeros(1,n_endor);
                HF_zx=zeros(1,n_endor);
                HF_zy=zeros(1,n_endor);

                for m=1:n_endor

                    HF_zz(m)=(sin(theta))^2*(cos(phi))^2*A(3*m-2,1)...
                        +(sin(theta))^2*(sin(phi))^2*A(3*m-1,2)...
                        +(cos(theta))^2*A(3*m,3)...
                        +2*(sin(theta))^2*sin(phi)*cos(phi)*A(3*m-2,2) ...
                        +2*sin(theta)*cos(theta)*cos(phi)*A(3*m-2,3)...
                        +2*sin(theta)*cos(theta)*sin(phi)*A(3*m-1,3);
                    %A11+A22+A33+A12+A13+A23

                    if parameters.Bterm==1
                        HF_zx(m)=(A(3*m-2,1)*sin(theta)*cos(phi)...
                            +A(3*m-1,1)*sin(theta)*sin(phi)...
                            +A(3*m,1)*cos(theta))*cos(theta)*cos(phi)...
                            +(A(3*m-2,2)*sin(theta)*cos(phi)...
                            +A(3*m-1,2)*sin(theta)*sin(phi)...
                            +A(3*m,2)*cos(theta))*cos(theta)*sin(phi)...
                            -(A(3*m-2,3)*sin(theta)*cos(phi)...
                            +A(3*m-1,3)*sin(theta)*sin(phi)...
                            +A(3*m,3)*cos(theta))*sin(theta);

                        HF_zy(m)=-A(3*m-2,1)*sin(theta)*cos(phi)*sin(phi)...
                            -A(3*m-1,1)*sin(theta)*sin(phi)^2 ...
                            -A(3*m,1)*cos(theta)*sin(phi)...
                            +A(3*m-2,2)*sin(theta)*cos(phi)^2 ...
                            +A(3*m-1,2)*sin(theta)*sin(phi)*cos(phi)...
                            +A(3*m,2)*cos(theta)*cos(phi);
                    end
                end

                % NQI ENDOR
                NQI_zz=zeros(1,n_endor);
                NQI=zeros(n_endor,3,3);

                for m=1:n_endor
                    if parameters.nqi_active==true

                        NQI_zz(m)=(sin(theta))^2*(cos(phi))^2*Q(3*m-2,1)...
                            +(sin(theta))^2*(sin(phi))^2*Q(3*m-1,2)...
                            +(cos(theta))^2*Q(3*m,3)...
                            +2*(sin(theta))^2*sin(phi)*cos(phi)*Q(3*m-2,2)...
                            +2*sin(theta)*cos(theta)*cos(phi)*Q(3*m-2,3)...
                            +2*sin(theta)*cos(theta)*sin(phi)*Q(3*m-1,3);
                        % Q11,Q22,Q33,Q12,Q13,Q23

                        % Rotation matrix into lab system
                        R1=zeros(3);
                        R1(1,1)=cos(theta)*cos(phi);
                        R1(1,2)=cos(theta)*sin(phi);
                        R1(1,3)=-sin(theta);
                        R1(2,1)=-sin(phi);
                        R1(2,2)=cos(phi);
                        R1(2,3)=0;
                        R1(3,1)=sin(theta)*cos(phi);
                        R1(3,2)=sin(theta)*sin(phi);
                        R1(3,3)=cos(theta);

                        Q_ENDOR=parameters.nqi_matrix;
                        qq2=Q_ENDOR((m-1)*3+1:(m-1)*3+3,:);
                        Y2=R1*qq2*R1';
                        NQI(m,:,:)=Y2;

                    else
                        NQI_zz(m)=0;
                        NQI (m,:,:)=0;
                    end
                end

                % CS ENDOR
                CS_zz=zeros(1,n_endor);
                if parameters.cs_active==true
                    for m=1:n_endor
                        CS_zz(m)=(sin(theta))^2*(cos(phi))^2*CS(3*m-2,1)...
                            +(sin(theta))^2*(sin(phi))^2*CS(3*m-1,2)...
                            +(cos(theta))^2*CS(3*m,3)...
                            +2*(sin(theta))^2*sin(phi)*cos(phi)*CS(3*m-2,2)...
                            +2*sin(theta)*cos(theta)*cos(phi)*CS(3*m-2,3)...
                            +2*sin(theta)*cos(theta)*sin(phi)*CS(3*m-1,3);
                    end
                end

                % dipolar SSC ENDOR
                D_zz=zeros(1,n_endor);

                if parameters.dipolar_active==true
                    D_zz=zeros(size(D,1)/3);
                    for m=1:size(D,1)/3
                        D_zz(m)=(sin(theta))^2*(cos(phi))^2*D(3*m-2,1)...
                            +(sin(theta))^2*(sin(phi))^2*D(3*m-1,2)...
                            +(cos(theta))^2*D(3*m,3)...
                            +2*(sin(theta))^2*sin(phi)*cos(phi)*D(3*m-2,2)...
                            +2*sin(theta)*cos(theta)*cos(phi)*D(3*m-2,3)...
                            +2*sin(theta)*cos(theta)*sin(phi)*D(3*m-1,3);
                    end
                end

                % HF EPR
                HF_zz_EPR=zeros(1,Ni_EPR);
                if Ni_EPR>0
                    for m=1:Ni_EPR
                        HF_zz_EPR(m)=(sin(theta))^2*(cos(phi))^2*A_EPR(3*m-2,1)...
                            +(sin(theta))^2*(sin(phi))^2*A_EPR(3*m-1,2)...
                            +(cos(theta))^2*A_EPR(3*m,3)...
                            +2*(sin(theta))^2*sin(phi)*cos(phi)*A_EPR(3*m-2,2)...
                            +2*sin(theta)*cos(theta)*cos(phi)*A_EPR(3*m-2,3)...
                            +2*sin(theta)*cos(theta)*sin(phi)*A_EPR(3*m-1,3);
                    end
                end

                fieldAxis=paramsEPR('fieldAxis');

                if Ni_EPR>0
                    I_EPR=parameters.epr_spin_numbers;
                    mI=kehl_ori_field_build_space(I_EPR);
                    for epr_idx=1:Ni_EPR
                        % EPR resonance with hyperfine coupling only
                        E=Beff+mI*(HF_zz_EPR(epr_idx)/constants("CONST1")*10^(-10))';
                    end
                    bin=(round((E-fieldAxis(1))/parameters.field_step_t)+1);
                else
                    E=Beff;
                    bin=(round((E-fieldAxis(1))/parameters.field_step_t)+1);
                end

                for p=1:length(bin)

                    % Add resonances to EPR spectrum
                    tmp_epr(bin(p))=tmp_epr(bin(p))+1;
                end

                offsets=kehl_offsets(constants,parameters,parameters.operator_spin_system,paramsENDOR,B,geff,HF_zz,NQI_zz);

                DeltaB=0;
                scalefactor=0;

                for root_idx=1:length(E)

                    % Magnetic field offset in T
                    DeltaB=fieldAxis(bin(root_idx))-parameters.field_t;

                    if isfield(parameters,'pulse_file')
                        scalefactor=kehl_ori_field_pulsescale(parameters,DeltaB,constants("CONST1"));
                    else

                        if abs(DeltaB)<=abs(parameters.nwidth*(W1/(constants("CONST1")*1e10)))

                            %exp('res_EN')
                            if abs(HF_zz(m))>1
                                Scale=(DeltaB^2-(W1/(constants("CONST1")*1e10))^2)/(DeltaB^2+(W1/(constants("CONST1")*1e10))^2);
                                Scale=(1-Scale)/2;
                            else
                                Scale=0;
                            end
                        else
                            Scale=0;
                        end
                        scalefactor=Scale;
                    end

                    %Select only those parameters, for which scalefactor > 0
                    if scalefactor>1e-3
                        or=or+1;
                        geff_sel(or)=geff(1);
                        B_sel(or)=Beff;
                        HF_zz_sel(or,:)=HF_zz(:);
                        HF_zy_sel(or,:)=HF_zy(:);
                        HF_zx_sel(or,:)=HF_zx(:);
                        NQI_zz_sel(or,:)=NQI_zz(:);
                        NQI_sel(or,:,:,:)=NQI(:,:,:);
                        CS_zz_sel(or,:)=CS_zz(:);
                        D_zz_sel(or,:)=D_zz(:);
                        euler_sel(or,:)=euler_angles;

                        S_sel(or)=scalefactor;
                        offsets_sel(or,:)=offsets(:);
                    end
                end
            end
            epr_amp=epr_amp+tmp_epr;
        end

        % only one orientation for single crystal calculation
    elseif parameters.powder==false
        or=or+1;
        geff=parameters.g_iso;

        %effective B field for given theta and phi
        Beff=parameters.mw_freq_hz*constants('H')/(constants('MU_B')*geff);

        HF_zz=zeros(1,n_endor);
        HF_zx=zeros(1,n_endor);
        HF_zy=zeros(1,n_endor);
        NQI_zz=zeros(1,n_endor);
        NQI=zeros(n_endor,3,3);

        CS_zz=zeros(1,n_endor);

        for m=1:n_endor

            HF_zz(m)=A(3*m,3);
            if parameters.Bterm==1
                HF_zy(m)=A(3*m,2);
                HF_zx(m)=A(3*m,1);
            end

            if parameters.nqi_active==true
                NQI_zz(m)=Q(3*m,3);
                NQI(m,:,:)=Q;
            end

            % CS ENDOR
            if parameters.cs_active==true
                CS_zz(m)=CS(3*m,3);
            end

        end

        % dipolar SSC ENDOR
        D_zz=zeros(1,n_endor);

        if parameters.dipolar_active==true
            D_zz=zeros(1,size(D,1)/3);
            for m=1:size(D,1)/3
                D_zz(m)=D(3*m,3);
            end
        end

        offsets_tmp=kehl_offsets(constants,parameters,parameters.operator_spin_system,paramsENDOR,Beff,geff,HF_zz,NQI_zz);

        sel_I=parameters.sel_I;
        if parameters.n_spin_systems>1
            s=size(offsets_tmp,2)/n_endor;
            for n=1:n_endor
                offsets(n)=offsets_tmp(sel_I+(n-1)*s);
            end
        else
            s=size(offsets_tmp,2)/n_endor;
            offsets=offsets_tmp(((sel_I-1)*s+1):(sel_I*s));
        end

        scalefactor=1;
        geff_sel(or)=geff(1);
        B_sel(or)=Beff;
        HF_zz_sel(or,:)=HF_zz(:);
        HF_zy_sel(or,:)=HF_zy(:);
        HF_zx_sel(or,:)=HF_zx(:);
        NQI_zz_sel(or,:)=NQI_zz(:);
        NQI_sel(or,:,:,:)=NQI(:,:,:);
        CS_zz_sel(or,:)=CS_zz(:);
        D_zz_sel(or,:)=D_zz(:);
        euler_sel(or,:)=[0 0 0];

        S_sel(or)=scalefactor;
        offsets_sel(or,:)=offsets(:);
    end
    % Build output map
    EPR=containers.Map;
    EPR("geff_sel")=geff_sel;
    EPR("B_sel")=B_sel;
    EPR("HF_zz_sel")=HF_zz_sel;
    EPR("HF_zy_sel")=HF_zy_sel;
    EPR("HF_zx_sel")=HF_zx_sel;
    EPR("CS_zz_sel")=CS_zz_sel;
    EPR("D_zz_sel")=D_zz_sel;
    EPR("euler_sel")=euler_sel;

    EPR("NQI_zz_sel")=NQI_zz_sel;
    EPR("NQI_sel")=NQI_sel;

    EPR("S_sel")=S_sel;
    EPR("EPR_amp")=epr_amp;

    EPR("offsets")=offsets_sel;
end
function grumble(spin_system,parameters)
    if (~isstruct(spin_system))||(~isfield(spin_system,'bas'))||(~isfield(spin_system,'comp'))
        error('spin_system must be a Spinach spin system structure.');
    end
    if ~isstruct(parameters)
        error('parameters must be a structure.');
    end
end

