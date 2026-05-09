% Frequency-domain EPR orientation selection for Kehl ENDOR. Syntax:
%
%      EPR=kehl_ori_freq(spin_system,parameters)
%
% Parameters:
%
%   spin_system      - Spinach spin system structure.
%   parameters       - Kehl ENDOR context parameter structure.
%
% Outputs:
%
%   EPR              - map containing selected frequency-domain orientations.
%
% February 2024 A. Kehl (akehl@gwdg.de)
% May 2026 Spinach integration
%
% <https://spindynamics.org/wiki/index.php?title=kehl_ori_freq.m>

function EPR=kehl_ori_freq(spin_system,parameters)

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
    Ni_EPR=parameters.n_epr;
    SzEPR=full(operator(spin_system,'Lz',parameters.electron_spin_idx));
    SxEPR=full(operator(spin_system,'Lx',parameters.electron_spin_idx));
    if parameters.epr_nuclei_active==true
        A_EPR=parameters.epr_hfc_matrix;
        Q_EPR=parameters.epr_nqi_matrix;
        g_N_EPR=parameters.epr_gamma_hz_t;
        IzEPR=cell(1,Ni_EPR);
        IxEPR=cell(1,Ni_EPR);
        IyEPR=cell(1,Ni_EPR);
        for n=1:Ni_EPR
            spin_idx=parameters.epr_spins(n);
            IzEPR{n}=full(operator(spin_system,'Lz',spin_idx));
            IxEPR{n}=full(operator(spin_system,'Lx',spin_idx));
            IyEPR{n}=full(operator(spin_system,'Ly',spin_idx));
        end
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
    offsets_sel=[];

    S_sel=[];

    nor=0;

    % Loop over powder orientations
    if parameters.powder==true
        if isfield(parameters,'excite_width')
            W1=parameters.excite_width;
        else
            W1=parameters.pulse_width*constants("CONST1")*1e10;
        end

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
                Beff=(parameters.field_t*parameters.g_iso)/geff(1);
                B=parameters.field_t;

                % Resonance frequency for geff at ObsField in GHz

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

                grot=R1*g*R1';

                % ENDOR values
                HF_zz=zeros(1,n_endor);
                HF_zx=zeros(1,n_endor);
                HF_zy=zeros(1,n_endor);
                NQI_zz=zeros(1,n_endor);
                NQI=zeros(n_endor,3,3);
                CS_zz=zeros(1,n_endor);
                D_zz=zeros(1,n_endor);

                for m=1:n_endor
                    % Hyperfine ENDOR projection
                    HF_zz(m)=(sin(theta))^2*(cos(phi))^2*A(3*m-2,1)...
                        +(sin(theta))^2*(sin(phi))^2*A(3*m-1,2)...
                        +(cos(theta))^2*A(3*m,3)...
                        +2*(sin(theta))^2*sin(phi)*cos(phi)*A(3*m-2,2) ...
                        +2*sin(theta)*cos(theta)*cos(phi)*A(3*m-2,3)...
                        +2*sin(theta)*cos(theta)*sin(phi)*A(3*m-1,3);

                    % Off-diagonal hyperfine projections for B-term correction
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

                    % NQI ENDOR
                    if parameters.nqi_active==true

                        NQI_zz(m)=(sin(theta))^2*(cos(phi))^2*Q(3*m-2,1)...
                            +(sin(theta))^2*(sin(phi))^2*Q(3*m-1,2)...
                            +(cos(theta))^2*Q(3*m,3)...
                            +2*(sin(theta))^2*sin(phi)*cos(phi)*Q(3*m-2,2)...
                            +2*sin(theta)*cos(theta)*cos(phi)*Q(3*m-2,3)...
                            +2*sin(theta)*cos(theta)*sin(phi)*Q(3*m-1,3);

                        % Full quadrupolar tensor in the laboratory frame
                        qq2=Q((m-1)*3+1:(m-1)*3+3,:);
                        Y2=R1*qq2*R1';
                        NQI(m,:,:)=Y2;
                    end

                    % CS ENDOR
                    if parameters.cs_active==true
                        CS_zz(m)=(sin(theta))^2*(cos(phi))^2*CS(3*m-2,1)...
                            +(sin(theta))^2*(sin(phi))^2*CS(3*m-1,2)...
                            +(cos(theta))^2*CS(3*m,3)...
                            +2*(sin(theta))^2*sin(phi)*cos(phi)*CS(3*m-2,2)...
                            +2*sin(theta)*cos(theta)*cos(phi)*CS(3*m-2,3)...
                            +2*sin(theta)*cos(theta)*sin(phi)*CS(3*m-1,3);
                    end

                end

                % Nuclear dipolar-dipolar coupling ENDOR projection
                if parameters.dipolar_active==true
                    D_zz=zeros(1,size(D,1)/3);
                    for m=1:size(D,1)/3
                        D_zz(m)=(sin(theta))^2*(cos(phi))^2*D(3*m-2,1)...
                            +(sin(theta))^2*(sin(phi))^2*D(3*m-1,2)...
                            +(cos(theta))^2*D(3*m,3)...
                            +2*(sin(theta))^2*sin(phi)*cos(phi)*D(3*m-2,2)...
                            +2*sin(theta)*cos(theta)*cos(phi)*D(3*m-2,3)...
                            +2*sin(theta)*cos(theta)*sin(phi)*D(3*m-1,3);
                    end
                end

                % EPR tensors
                HF_EPR=zeros(3,3,Ni_EPR);
                NQI_EPR=zeros(3,3,Ni_EPR);

                if Ni_EPR>0
                    for m=1:Ni_EPR
                        % Hyperfine EPR tensor
                        hf2=A_EPR((m-1)*3+1:(m-1)*3+3,:);
                        X2=R1*hf2*R1';
                        HF_EPR(:,:,m)=X2;

                        % Quadrupolar EPR tensor
                        if parameters.epr_nqi_active==true
                            qq2=Q_EPR((m-1)*3+1:(m-1)*3+3,:);
                            Y2=R1*qq2*R1';
                            NQI_EPR(:,:,m)=Y2;
                        end

                    end
                end

                H_EZ_EPR=constants("MU_B")*grot(3,3)*B/constants("H")*SzEPR;
                H_NZ_EPR=zeros(size(SzEPR));
                H_HF_EPR=zeros(size(SzEPR));
                H_NQI_EPR=zeros(size(SzEPR));

                if parameters.epr_nuclei_active>0
                    for mm=1:Ni_EPR
                        H_NZ_EPR=H_NZ_EPR+g_N_EPR(mm)*IzEPR{mm};
                        H_HF_EPR=H_HF_EPR+IzEPR{mm}*HF_EPR(3,3,mm)*SzEPR+IxEPR{mm}*(HF_EPR(1,3,mm)^2+HF_EPR(2,3,mm)^2)^0.5*SzEPR;
                        H_NQI_EPR=H_NQI_EPR+NQI_EPR(1,1,mm)*IxEPR{mm}*IxEPR{mm};
                        H_NQI_EPR=H_NQI_EPR+NQI_EPR(1,2,mm)*IxEPR{mm}*IyEPR{mm};
                        H_NQI_EPR=H_NQI_EPR+NQI_EPR(1,3,mm)*IxEPR{mm}*IzEPR{mm};
                        H_NQI_EPR=H_NQI_EPR+NQI_EPR(2,1,mm)*IyEPR{mm}*IxEPR{mm};
                        H_NQI_EPR=H_NQI_EPR+NQI_EPR(2,2,mm)*IyEPR{mm}*IyEPR{mm};
                        H_NQI_EPR=H_NQI_EPR+NQI_EPR(2,3,mm)*IyEPR{mm}*IzEPR{mm};
                        H_NQI_EPR=H_NQI_EPR+NQI_EPR(3,1,mm)*IzEPR{mm}*IxEPR{mm};
                        H_NQI_EPR=H_NQI_EPR+NQI_EPR(3,2,mm)*IzEPR{mm}*IyEPR{mm};
                        H_NQI_EPR=H_NQI_EPR+NQI_EPR(3,3,mm)*IzEPR{mm}*IzEPR{mm};
                    end
                end

                H_S_EPR=H_EZ_EPR-H_NZ_EPR+H_HF_EPR+H_NQI_EPR;

                % Diagonalize and calculate eigenvalues of the EPR Hamiltonian

                % E ... Hamiltonian Elements, V... eigenvectors
                [V_EPR,E_EPR]=eig(H_S_EPR);

                % round necessary for following calculation
                V_EPR=round((V_EPR),9);
                E_EPR=real(diag(E_EPR));

                % Calculate transitions between all elements
                trans_prob_EPR=[];
                freq_EPR=[];
                q=1;
                for x=1:length(V_EPR)
                    for y=x+1:length(V_EPR)
                        % Electron transition probability in the eigenbasis
                        trans_prob_EPR(q)=abs(round((V_EPR(:,x))'*SxEPR*(V_EPR(:,y)),9))^2;

                        % Transition frequency in Hz
                        freq_EPR(q)=abs(E_EPR(x)-E_EPR(y));

                        % Select EPR transitions with frequency threshold value (1 GHz)
                        if (freq_EPR(q)<1e9)||(trans_prob_EPR(q)<0.1)
                            trans_prob_EPR(q)=0;
                        end
                        q=q+1;
                    end
                end

                freq_EPR=freq_EPR(logical(trans_prob_EPR))   ;
                trans_prob_EPR=trans_prob_EPR(logical(trans_prob_EPR));

                % EPR Frequency Spectrum
                for p=1:length(freq_EPR)

                    % Scale resonance to frequency-axis bin
                    bin_freq=round((freq_EPR-parameters.epr_freq_min_hz)/parameters.epr_freq_step_hz)+1;
                    if trans_prob_EPR(p)>0
                        tmp_epr(bin_freq(p))=tmp_epr(bin_freq(p))+trans_prob_EPR(p);
                        DeltaOm=freq_EPR(p)-parameters.mw_freq_hz;

                        if isfield(parameters,'pulse_file')
                            SF=kehl_ori_freq_pulsescale(parameters,DeltaOm)/2;
                        else
                            SF=((DeltaOm^2)-(W1^2))/((DeltaOm^2)+(W1^2));
                            SF=(1-SF)/2;
                        end

                        if SF>1e-3

                            % including transition probability!
                            scalefactor=SF;
                        else
                            scalefactor=0;
                        end
                    end

                    %Select only those parameters, for which scalefactor > 0

                    offsets=kehl_offsets(constants,parameters,spin_system,paramsENDOR,B,geff,HF_zz,NQI_zz);

                    if (scalefactor>0)&&(trans_prob_EPR(p)>0)

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

        % single orientation for single crystal
    elseif parameters.powder==false
        or=or+1;
        geff=parameters.g_iso;

        %effective B field for given theta and phi
        B=(parameters.field_t*parameters.g_iso)/geff(1);

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

        % nuclear dipolar-dipolar coupling ENDOR
        D_zz=zeros(1,n_endor);
        if parameters.dipolar_active==true
            D_zz=zeros(1,size(D,1)/3);
            for m=1:size(D,1)/3
                D_zz(m)=D(3*m,3);
            end
        end

        offsets_tmp=kehl_offsets(constants,parameters,spin_system,paramsENDOR,B,geff,HF_zz,NQI_zz);

        sel_I=parameters.sel_I;
        s=size(offsets_tmp,2)/n_endor;
        offsets=offsets_tmp(((sel_I-1)*s+1):(sel_I*s));

        scalefactor=1;
        geff_sel(or)=geff(1);
        B_sel(or)=B;
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

    % build output Map
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

