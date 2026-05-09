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
    paramsEPR=parameters.paramsEPR;
    Ntheta=parameters.Nang;
    Nphimax=parameters.Nang;
    g=parameters.g_matrix;
    A=parameters.hfc_matrix;
    Q=parameters.nqi_matrix;
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

    % Loop over powder orientations
    if parameters.powder==true
        W1=parameters.excite_width;

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
                Beff=(parameters.static_field*parameters.g_iso)/geff(1);
                B=parameters.static_field;

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

                % Select microwave-active EPR transitions in Liouville space
                [freq_EPR,trans_prob_EPR]=kehl_epr_transitions(...
                    spin_system,parameters,euler_angles,'frequency',B);
                hit_list=(freq_EPR<parameters.epr_freq_min)|(trans_prob_EPR<0.1);
                freq_EPR(hit_list)=[];
                trans_prob_EPR(hit_list)=[];

                % EPR Frequency Spectrum
                for p=1:length(freq_EPR)

                    % Scale resonance to frequency-axis bin
                    scalefactor=0;
                    bin_freq=round((freq_EPR(p)-parameters.epr_freq_min)/...
                        parameters.epr_freq_step)+1;
                    if (bin_freq>=1)&&(bin_freq<=paramsEPR("Npts"))&&...
                            (trans_prob_EPR(p)>0)
                        tmp_epr(bin_freq)=tmp_epr(bin_freq)+trans_prob_EPR(p);
                        DeltaOm=freq_EPR(p)-parameters.mw_freq;

                        if isfield(parameters,'pulse_file')
                            SF=kehl_ori_freq_pulsescale(parameters,DeltaOm)/2;
                        else
                            SF=((DeltaOm^2)-(W1^2))/((DeltaOm^2)+(W1^2));
                            SF=(1-SF)/2;
                        end

                        if SF>1e-3

                            % Include transition probability in the EPR spectrum
                            scalefactor=SF;
                        end
                    end

                    % Select only transitions with non-negligible excitation
                    offsets=kehl_offsets(parameters,spin_system,B,euler_angles);

                    if (scalefactor>0)&&(trans_prob_EPR(p)>0)

                        or=or+1;
                        geff_sel(or)=geff(1); %#ok<AGROW>
                        B_sel(or)=Beff; %#ok<AGROW>
                        HF_zz_sel(or,:)=HF_zz(:); %#ok<AGROW>
                        HF_zy_sel(or,:)=HF_zy(:); %#ok<AGROW>
                        HF_zx_sel(or,:)=HF_zx(:); %#ok<AGROW>
                        NQI_zz_sel(or,:)=NQI_zz(:); %#ok<AGROW>
                        NQI_sel(or,:,:,:)=NQI(:,:,:); %#ok<AGROW>
                        CS_zz_sel(or,:)=CS_zz(:); %#ok<AGROW>
                        D_zz_sel(or,:)=D_zz(:); %#ok<AGROW>
                        euler_sel(or,:)=euler_angles; %#ok<AGROW>
                        S_sel(or)=scalefactor; %#ok<AGROW>
                        offsets_sel(or,:)=offsets(:); %#ok<AGROW>
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
        B=(parameters.static_field*parameters.g_iso)/geff(1);

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

        offsets=kehl_offsets(parameters,spin_system,B,[0 0 0]);

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
    if ~strcmp(spin_system.bas.formalism,'zeeman-liouv')
        error('spin_system must use zeeman-liouv formalism.');
    end
    if ~isstruct(parameters)
        error('parameters must be a structure.');
    end
end

