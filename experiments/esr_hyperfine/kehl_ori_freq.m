% calculates the selected EPR orientations in the frequency domain and the
% corresponding effective spin parameter values
%
% input parameters:
% constants: the Map containing the constants
% spinSys: the Map describing the spin system
% spinOps: the Map containing the spin operators
% paramsEPR: the Map containing the EPR parameters
% opt: the Map containing the optional parameters
% expt: the Map containing the experimental parameters
%
% output parameters:
% EPR: Map containing the information from the EPR experiment
%
% February 2024 A. Kehl (akehl@gwdg.de)
%

% get parameters from Maps

function EPR=kehl_ori_freq(constants,spinSys,spinOps,paramsEPR,paramsENDOR,opt,expt)

    % Check consistency
    grumble(constants,spinSys,spinOps,paramsEPR,paramsENDOR,opt,expt);
    Ntheta=opt("Nang");
    Nphimax=opt("Nang");
    g=spinSys("g");
    A=spinSys("A");
    Q=spinSys("Q");
    if spinSys("EPR_Nucs_used")==true
        Ni_EPR=spinSys("Ni_EPR");
        A_EPR=spinSys("A_EPR");
        Q_EPR=spinSys("Q_EPR");
        g_N_EPR=spinSys("g_N_EPR");
        I_EPR=spinSys("I_EPR");
        spinOpsEPR=kehl_spin_ops(spinSys("S"),I_EPR,Ni_EPR,Ni_EPR);
        SzEPR=spinOpsEPR("Sz");
        SxEPR=spinOpsEPR("Sx");
        IzEPR=spinOpsEPR("Iz");
        IxEPR=spinOpsEPR("Ix");
        IyEPR=spinOpsEPR("Iy");
    else
        Ni_EPR=0;
    end
    if spinSys("CS_used")==true
        CS=spinSys("CS");
    end
    if spinSys("D_used")==true
        D=spinSys("D");
    end
    Ni_ENDOR=spinSys("Ni_ENDOR");
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
    offsets_sel=[];

    S_sel=[];

    nor=0;

    % loop over orientations for powder pattern
    if opt("powder")==true
        if isKey(expt,'exciteWidth')
            W1=expt('exciteWidth');
        else
            W1=expt("pulsewidth")*constants("CONST1")*1e10;
        end

        for ii=1:Ntheta
            theta=ii*pi/Ntheta;
            Nphi=round(sin(theta)*Nphimax)*1;
            for jj=1:Nphi
                phi=(jj-1)*pi*2/(Nphi);

                % direction cosine vector
                dc=[cos(phi)*sin(theta) sin(phi)*sin(theta) cos(theta)];
                % Effective g-value for a given theta, and phi combination
                geff=(dc*g2*dc')^.5;

                % effective B field for given theta and phi
                Beff=(expt("Field")*spinSys("g_iso"))/geff(1);
                B=expt("Field");

                % Resonance frequency for geff at ObsField in GHz
                veff=geff*expt("Field")*9.27401*10^-24/(6.62607*10^-34);

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
                HF_zz=zeros(1,Ni_ENDOR);
                HF_zx=zeros(1,Ni_ENDOR);
                HF_zy=zeros(1,Ni_ENDOR);
                NQI_zz=zeros(1,Ni_ENDOR);
                NQI=zeros(Ni_ENDOR,3,3);
                CS_zz=zeros(1,Ni_ENDOR);
                D_zz=zeros(1,Ni_ENDOR);

                for m=1:Ni_ENDOR
                     % HF ENDOR
                     HF_zz(m)=(sin(theta))^2*(cos(phi))^2*A(3*m-2,1)...
+(sin(theta))^2*(sin(phi))^2*A(3*m-1,2)...
+(cos(theta))^2*A(3*m,3)...
+2*(sin(theta))^2*sin(phi)*cos(phi)*A(3*m-2,2) ...
+2*sin(theta)*cos(theta)*cos(phi)*A(3*m-2,3)...
+2*sin(theta)*cos(theta)*sin(phi)*A(3*m-1,3);
                        %A11+A22+A33+A12+A13+A23
                     if opt("Bterm")==1
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
                     if spinSys("Q_used")==true

                        NQI_zz(m)=(sin(theta))^2*(cos(phi))^2*Q(3*m-2,1)...
+(sin(theta))^2*(sin(phi))^2*Q(3*m-1,2)...
+(cos(theta))^2*Q(3*m,3)...
+2*(sin(theta))^2*sin(phi)*cos(phi)*Q(3*m-2,2)...
+2*sin(theta)*cos(theta)*cos(phi)*Q(3*m-2,3)...
+2*sin(theta)*cos(theta)*sin(phi)*Q(3*m-1,3);
                        % Q11,Q22,Q33,Q12,Q13,Q23

                        % Q_ENDOR = spinSys("Q");
                        qq2=Q((m-1)*3+1:(m-1)*3+3,:);
                        Y2=R1*qq2*R1';
                        NQI(m,:,:)=Y2;
                     end

                     % CS ENDOR
                     if spinSys("CS_used")==true
                        CS_zz(m)=(sin(theta))^2*(cos(phi))^2*CS(3*m-2,1)...
+(sin(theta))^2*(sin(phi))^2*CS(3*m-1,2)...
+(cos(theta))^2*CS(3*m,3)...
+2*(sin(theta))^2*sin(phi)*cos(phi)*CS(3*m-2,2)...
+2*sin(theta)*cos(theta)*cos(phi)*CS(3*m-2,3)...
+2*sin(theta)*cos(theta)*sin(phi)*CS(3*m-1,3);
                     end

                end
                % nulcear dipolar-dipolar coupling ENDOR
                 if spinSys("D_used")==true
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

                % EPR vlaues

                HF_EPR=zeros(3,3,Ni_EPR);
                NQI_EPR=zeros(3,3,Ni_EPR);
                % CS_EPR = zeros(3,3,Ni_EPR);

                if Ni_EPR>0
                    for m=1:Ni_EPR
                        % HF EPR
                        hf2=A_EPR((m-1)*3+1:(m-1)*3+3,:);
                        X2=R1*hf2*R1';
                        HF_EPR(:,:,m)=X2;


                         % Q EPR
                        if spinSys("EPR_Q_used")==true
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

               if spinSys('EPR_Nucs_used')>0
                   for mm=1:Ni_EPR
                        H_NZ_EPR=H_NZ_EPR+g_N_EPR(mm)*IzEPR{mm};
                        H_HF_EPR=H_HF_EPR+IzEPR{mm}*HF_EPR(3,3,mm)*SzEPR+IxEPR{mm}*(HF_EPR(1,3,mm)^2+HF_EPR(2,3,mm)^2)^0.5*SzEPR;
                        H_NQI_EPR=H_NQI_EPR+NQI_EPR(1,1,m)*IxEPR{mm}*IxEPR{mm};
                        H_NQI_EPR=H_NQI_EPR+NQI_EPR(1,2,m)*IxEPR{mm}*IyEPR{mm};
                        H_NQI_EPR=H_NQI_EPR+NQI_EPR(1,3,m)*IxEPR{mm}*IzEPR{mm};
                        H_NQI_EPR=H_NQI_EPR+NQI_EPR(2,1,m)*IyEPR{mm}*IxEPR{mm};
                        H_NQI_EPR=H_NQI_EPR+NQI_EPR(2,2,m)*IyEPR{mm}*IyEPR{mm};
                        H_NQI_EPR=H_NQI_EPR+NQI_EPR(2,3,m)*IyEPR{mm}*IzEPR{mm};
                        H_NQI_EPR=H_NQI_EPR+NQI_EPR(3,1,m)*IzEPR{mm}*IxEPR{mm};
                        H_NQI_EPR=H_NQI_EPR+NQI_EPR(3,2,m)*IzEPR{mm}*IyEPR{mm};
                        H_NQI_EPR=H_NQI_EPR+NQI_EPR(3,3,m)*IzEPR{mm}*IzEPR{mm};
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
                q=1;
                for x=1:length(V_EPR)
                    for y=x+1:length(V_EPR)
                        % Electron transition probability in the eigenbasis
                        trans_prob_EPR(q)=abs(round((V_EPR(:,x))'*SxEPR*(V_EPR(:,y)),9))^2;

                        %in Hz, Transition Frequency
                        freq_EPR(q)=abs(E_EPR(x)-E_EPR(y));

                        %Select EPR Transitions with frequency threshhold value (1 GHz)
                        if (freq_EPR(q)<1e9) || (trans_prob_EPR(q)<0.1)
                            trans_prob_EPR(q)=0;
                        end
                        q=q+1;
                    end
                end

                freq_EPR=freq_EPR(logical(trans_prob_EPR))   ;
                trans_prob_EPR=trans_prob_EPR(logical(trans_prob_EPR));

                % EPR_Resonances =[freq_EPR' trans_prob_EPR' grot(3,3)*ones(length(freq_EPR),1)];

                % EPR Frequency Spectrum
                for p=1:length(freq_EPR)

                    %scale resonance to freq axis bin
                    bin_freq=round((freq_EPR-expt("FreqMin"))/expt("FreqSteps"))+1;
                    if trans_prob_EPR(p)>0
                        tmp_epr(bin_freq(p))=tmp_epr(bin_freq(p))+trans_prob_EPR(p);
                        DeltaOm=freq_EPR(p)-expt("FreqMeas");

                        if isKey(expt,"pulse")
                            SF=pulsescale(expt,DeltaOm)/2;
                        else
                            SF=((DeltaOm^2)-(W1^2))/((DeltaOm^2)+(W1^2));
                            SF=(1-SF)/2;
                            nor=nor+1;
                        end

                        if SF>1e-3

                            % including transition probability!
                            scalefactor=SF;
                        else
                            scalefactor=0;
                        end
                    end


                    %Select only those parameters, for which scalefactor > 0

                    offsets=kehl_offsets(constants,spinSys,spinOps,paramsENDOR,B,geff,HF_zz,NQI_zz);


                    if (scalefactor>0) && (trans_prob_EPR(p)>0)

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
                        S_sel(or)=scalefactor;
                        offsets_sel(or,:)=offsets(:);
                    end
                end
            end
            epr_amp=epr_amp+tmp_epr;
        end

    % single orientation for single crystal
    elseif opt("powder")==false
        or=or+1;
        geff=spinSys("g_iso");

        %effective B field for given theta and phi
        B=(expt("Field")*spinSys("g_iso"))/geff(1);

        HF_zz=zeros(1,Ni_ENDOR);
        HF_zx=zeros(1,Ni_ENDOR);
        HF_zy=zeros(1,Ni_ENDOR);
        NQI_zz=zeros(1,Ni_ENDOR);
        NQI=zeros(Ni_ENDOR,3,3);


        for m=1:Ni_ENDOR

            HF_zz(m)=A(3*m,3);
            if opt("Bterm")==1
                HF_zy(m)=A(3*m,2);
                HF_zx(m)=A(3*m,1);
            end

            if spinSys("Q_used")==true
                NQI_zz(m)=Q(3*m,3);
                NQI(m,:,:)=Q;
            end

            % CS ENDOR
            CS_zz=zeros(1,Ni_ENDOR);
            if spinSys("CS_used")==true
                CS_zz(m)=CS(3*m,3);
            end
        end

        % nuclear dipolar-dipolar coupling ENDOR
        D_zz=zeros(1,Ni_ENDOR);
        if spinSys("D_used")==true
             D_zz=zeros(1,size(D,1)/3);
             for m=1:size(D,1)/3
                 D_zz(m)=D(3*m,3);
             end
         end

        offsets_tmp=kehl_offsets(constants,spinSys,spinOps,paramsENDOR,B,geff,HF_zz,NQI_zz);

        sel_I=opt("sel_I");
        if spinSys('N_SpinSys')>1
            s=size(offsets_tmp,2)/Ni_ENDOR;
            for n=1:Ni_ENDOR
                offsets(n)=offsets_tmp(sel_I+(n-1)*s);
            end
        else
            s=size(offsets_tmp,2)/Ni_ENDOR;
            offsets=offsets_tmp(((sel_I-1)*s+1):(sel_I*s));
        end

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

    EPR("NQI_zz_sel")=NQI_zz_sel;
    EPR("NQI_sel")=NQI_sel;

    EPR("S_sel")=S_sel;
    EPR("EPR_amp")=epr_amp;

    EPR("offsets")=offsets_sel;
end

function scalefactor=pulsescale(expt,DeltaOm)
    % calculates the scalefactor in dependence of the mw pulse's excitation
    % profile
    %
    % input parameters:
    % expt: the Map containing the experimental parameters
    % DeltaOm: offset of the actual freq from the resonance freq
    %
    % output parameters:
    % scalefactor: scalefactor for the specific orientation
    %
    % February 2024 A. Kehl (akehl@gwdg.de)
    %
    data=expt("pulse");
    %%
    pulse=fopen(data);

    pulse_data=textscan(pulse,'%f %f %f %f');
    fclose('all');
    pulse_x=pulse_data{2};
    dx=1000/(pulse_x(2)-pulse_x(1));

    n=1000;


    pulse_y=pulse_data{3}+1i*pulse_data{4};

    Y=abs(fftshift(fft(pulse_y(1:70),(n+1))));
    y=Y.^3;
    y=y/max(y);
    Y=Y/max(Y);
    if isKey(expt,"3pulses") && expt("3pulses")==true
       Y=y;
    end
%%

    %Delta omega in MHz
    DeltaOM=DeltaOm*2*pi*1e-6;
    binDOM=round(DeltaOM)+(dx/2+1);
    if binDOM>0 && binDOM<(n+1)
        scalefactor=Y(binDOM);
    else
        scalefactor=0;
    end
end

function grumble(constants,spinSys,spinOps,paramsEPR,paramsENDOR,opt,expt)
if ~isa(constants,'containers.Map')
    error('constants must be a containers.Map object.');
end
if (~isempty(spinSys))&&(~isa(spinSys,'containers.Map'))
    error('spinSys must be empty, or a containers.Map object.');
end
if ~isa(spinOps,'containers.Map')
    error('spinOps must be a containers.Map object.');
end
if ~isa(paramsEPR,'containers.Map')
    error('paramsEPR must be a containers.Map object.');
end
if ~isa(paramsENDOR,'containers.Map')
    error('paramsENDOR must be a containers.Map object.');
end
if ~isa(opt,'containers.Map')
    error('opt must be a containers.Map object.');
end
if ~isa(expt,'containers.Map')
    error('expt must be a containers.Map object.');
end
end

