% calculates the selected EPR orientations in the field domain and the
% corresponding effective spin parameter values
%
% input parameters:
% constants: the Map containing the constants
% spinSys: the Map describing the spin system
% spinOps: the Map containing the spin operators
% paramsEPR: the Map containing the EPR parameters
% parameters: structure containing simulation parameters
% expt: the Map containing the experimental parameters
%
% output parameters:
% EPR: Map containing the information from the EPR experiment
%
% February 2024 A. Kehl (akehl@gwdg.de)
%

% get parameters from Maps

function EPR=kehl_ori_field(constants,spinSys,spinOps,paramsEPR,paramsENDOR,parameters,expt)

    % Check consistency
    grumble(constants,spinSys,spinOps,paramsEPR,paramsENDOR,parameters,expt);
    Ntheta=parameters.Nang;
    Nphimax=parameters.Nang;
    g=spinSys("g");
    A=spinSys("A");
    Q=spinSys("Q");
    if spinSys("EPR_Nucs_used")==true
        Ni_EPR=spinSys("Ni_EPR");
        A_EPR=spinSys("A_EPR");
        Q_EPR=spinSys("Q_EPR");
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
    S_sel=[];
    offsets_sel=[];


    % for powder pattern
    if parameters.powder==true
        if isKey(expt,'exciteWidth')
            W1=expt('exciteWidth');
        else
            W1=expt("pulsewidth")*constants("CONST1")*1e10;
        end
        nor=0;
    % loop through orientations
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
            B=(expt("Field")*spinSys("g_iso"))/geff(1);
            Beff=expt('FreqMeas')*constants('H')/(constants('MU_B')*geff);

            % HF ENDOR
            HF_zz=zeros(1,Ni_ENDOR);
            HF_zx=zeros(1,Ni_ENDOR);
            HF_zy=zeros(1,Ni_ENDOR);

            for m=1:Ni_ENDOR

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
            NQI_zz=zeros(1,Ni_ENDOR);
            NQI=zeros(Ni_ENDOR,3,3);

            for m=1:Ni_ENDOR
                if spinSys("Q_used")==true

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

                    Q_ENDOR=spinSys("Q");
                    qq2=Q_ENDOR((m-1)*3+1:(m-1)*3+3,:);
                    Y2=R1*qq2*R1';
                    NQI(m,:,:)=Y2;

                else
                    NQI_zz(m)=0;
                    NQI (m,:,:)=0;
                end
            end

            % CS ENDOR
            CS_zz=zeros(1,Ni_ENDOR);
            if spinSys("CS_used")==true
                for m=1:Ni_ENDOR
                    CS_zz(m)=(sin(theta))^2*(cos(phi))^2*CS(3*m-2,1)...
+(sin(theta))^2*(sin(phi))^2*CS(3*m-1,2)...
+(cos(theta))^2*CS(3*m,3)...
+2*(sin(theta))^2*sin(phi)*cos(phi)*CS(3*m-2,2)...
+2*sin(theta)*cos(theta)*cos(phi)*CS(3*m-2,3)...
+2*sin(theta)*cos(theta)*sin(phi)*CS(3*m-1,3);
                end
            end

            % dipolar SSC ENDOR
            D_zz=zeros(1,Ni_ENDOR);

             if spinSys("D_used")==true
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
                I_EPR=spinSys("I_EPR");
                mI=BuildSpace(I_EPR);
                for i=1:Ni_EPR
                        % EPR resonance with hyperfine coupling only
                        E=Beff+mI*(HF_zz_EPR(i)/constants("CONST1")*10^(-10))';
                end
                bin=(round((E-fieldAxis(1))/expt("deltaField"))+1);
            else
                E=Beff;
                bin=(round((E-fieldAxis(1))/expt("deltaField"))+1);
            end

            for p=1:length(bin)

                % Add resonances to EPR spectrum
                tmp_epr(bin(p))=tmp_epr(bin(p))+1;
            end

            offsets=kehl_offsets(constants,spinSys,spinOps,paramsENDOR,B,geff,HF_zz,NQI_zz);

            DeltaB=0;
            scalefactor=0;

            for l=1:length(E)
                nor=nor+1;

                % Magnetic field offset in T
                DeltaB=fieldAxis(bin(l))-expt("Field");

                if isKey(expt,"pulse")
                    scalefactor=pulsescale(expt,DeltaB,constants("CONST1"));
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
        geff=spinSys("g_iso");

        %effective B field for given theta and phi
        Beff=expt('FreqMeas')*constants('H')/(constants('MU_B')*geff);

        HF_zz=zeros(1,Ni_ENDOR);
        HF_zx=zeros(1,Ni_ENDOR);
        HF_zy=zeros(1,Ni_ENDOR);
        NQI_zz=zeros(1,Ni_ENDOR);
        NQI=zeros(Ni_ENDOR,3,3);


        CS_zz=zeros(1,Ni_ENDOR);

        for m=1:Ni_ENDOR

            HF_zz(m)=A(3*m,3);
            if parameters.Bterm==1
                HF_zy(m)=A(3*m,2);
                HF_zx(m)=A(3*m,1);
            end

            if spinSys("Q_used")==true
                NQI_zz(m)=Q(3*m,3);
                NQI(m,:,:)=Q;
            end

            % CS ENDOR
            if spinSys("CS_used")==true
                CS_zz(m)=CS(3*m,3);
            end

        end

        % dipolar SSC ENDOR
        D_zz=zeros(1,Ni_ENDOR);

        if spinSys("D_used")==true
             D_zz=zeros(1,size(D,1)/3);
             for m=1:size(D,1)/3
                 D_zz(m)=D(3*m,3);
             end
         end


        offsets_tmp=kehl_offsets(constants,spinSys,spinOps,paramsENDOR,Beff,geff,HF_zz,NQI_zz);

        sel_I=parameters.sel_I;
        if spinSys('N_SpinSys')>1
            s=size(offsets_tmp,2)/Ni_ENDOR;
            for n=1:Ni_ENDOR
                offsets(n)=offsets_tmp(sel_I+(n-1)*s);
            end
        else
            size(offsets_tmp,2)
            s=size(offsets_tmp,2)/Ni_ENDOR;
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

function scalefactor=pulsescale(expt,DeltaB,CONST1)
    % calculates the scalefactor in dependence of the mw pulse's excitation
    % profile
    %
    % input parameters:
    % expt: the Map containing the experimental parameters
    % DeltaB: offset of the actual field from the effective resonance field
    % CONST1: constant to convert between field and frequency
    %
    % output parameters:
    % scalefactor: scalefactor for the specific orientation
    %
    % February 2024 A. Kehl (akehl@gwdg.de)
    %

    data=expt("pulse");
    pulse=fopen(data);

    pulse_data=textscan(pulse,'%f %f %f %f');
    fclose('all');
    pulse_x=pulse_data{2};
    dx=1000/(pulse_x(2)-pulse_x(1));

    n=1000;


    pulse_y=pulse_data{3}+1i*pulse_data{4};

    Y=abs(fftshift(fft(pulse_y,(n+1))));
    y=Y.^3;
    y=y/max(y);
    Y=Y/max(Y);
    if isKey(expt,"3pulses") && expt("3pulses")==true
       Y=y;
    end

    %Delta omega/2pi in MHz
    DeltaOM=DeltaB*CONST1*10^4*2*pi;

    binDOM=round(DeltaOM)+(dx/2+1);
    if binDOM>0 && binDOM<(n+1)
        scalefactor=Y(binDOM);
    else
        scalefactor=0;
    end

end

function M=BuildSpace(S)
    % builds up the space of basis vectors in columns,
    % script by M. Bennati
    %
    % input parameters:
    % S: spin quantum number
    %
    % output parameters:
    % M: basis vectors in columns
    %
    % February 2024 A. Kehl (akehl@gwdg.de)
    %

    % Initialize the space to a null
    M=[];

    % Initialize the dimension to one
    dim=1;

    % Loop over each spin
    for j=1:length(S)

      % Initialize the holder
      temp=[];

      % Loop over the z components
      for k=-S(j):S(j)

        %   adding block of the old space
        temp=[temp; M k*ones(dim,1)];

      %   plus a column of the z component
      end

      % Update M
      M=temp;

      % The dimensionality is the number
      dim=size(M,1);
                                            %   of rows in M
    end
end

function grumble(constants,spinSys,spinOps,paramsEPR,paramsENDOR,parameters,expt)
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
if ~isstruct(parameters)
    error('parameters must be a structure.');
end
if ~isa(expt,'containers.Map')
    error('expt must be a containers.Map object.');
end
end

