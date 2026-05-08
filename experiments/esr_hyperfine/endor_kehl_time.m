%ENDOR_KEHL_TIME time-domain Mims ENDOR pulse sequence for the Kehl ENDOR context.
%
%   ENDOR_AMP=ENDOR_KEHL_TIME(SPIN_SYSTEM,PARAMETERS,H,R,K) is called by
%   endor_kehl_context.m. H, R, and K are present for Spinach experiment
%   signature compatibility; this sequence uses orientation-selected
%   effective Hamiltonian data prepared by the context.


function endor_amp=endor_kehl_time(spin_system,parameters,H,R,K)

% Check consistency
grumble(spin_system,parameters,H,R,K);
if parameters.Relax==true
    endor_amp=kehl_time_rlx(spin_system,parameters);
else
    endor_amp=kehl_time_calc(spin_system,parameters);
end
end

% Consistency enforcement
function grumble(spin_system,parameters,H,R,K)
if (~isstruct(spin_system))||(~isfield(spin_system,'bas'))||(~isfield(spin_system,'comp'))
    error('spin_system must be a Spinach spin system structure.');
end
if ~isstruct(parameters)
    error('parameters must be a structure.');
end
if (~isempty(H))&&(~isnumeric(H))
    error('H must be empty or numeric.');
end
if (~isempty(R))&&(~isnumeric(R))
    error('R must be empty or numeric.');
end
if (~isempty(K))&&(~isnumeric(K))
    error('K must be empty or numeric.');
end
end

function endor_amp=kehl_time_rlx(spin_system,parameters)

    % Check consistency
    kehl_time_rlx_grumble(spin_system,parameters);

    % Unpack context data
    constants=parameters.constants;
    expt=parameters.expt;
    paramsENDOR=parameters.paramsENDOR;
    EPR=parameters.epr;
    operator_spin_system=parameters.operator_spin_system;
    n_endor=parameters.n_endor;
    n_spin_systems=parameters.n_spin_systems;
    I=parameters.endor_spin_numbers;

    % Obtain electron operators directly from Spinach
    Sx=full(operator(operator_spin_system,'Lx',1));
    Sy=full(operator(operator_spin_system,'Ly',1));
    Sz=full(operator(operator_spin_system,'Lz',1));

    % Obtain nuclear operators directly from Spinach
    Ix=cell(1,operator_spin_system.comp.nspins-1);
    Iy=cell(1,operator_spin_system.comp.nspins-1);
    Iz=cell(1,operator_spin_system.comp.nspins-1);
    for n=1:numel(Ix)
        spin_idx=n+1;
        Ix{n}=full(operator(operator_spin_system,'Lx',spin_idx));
        Iy{n}=full(operator(operator_spin_system,'Ly',spin_idx));
        Iz{n}=full(operator(operator_spin_system,'Lz',spin_idx));
    end
    [Sx_D,Sy_D,Sz_D,Ix_D,Iy_D,Iz_D]=kehl_diag_ops(Sx,Sy,Sz,Ix,Iy,Iz,n_endor);



    t=expt("t");
    Nint=8;

    if n_spin_systems>1
       nucs=n_endor;
        for i=1:n_endor-1
            Ix{i+1}=Ix{1};
            Iy{i+1}=Iy{1};
            Iz{i+1}=Iz{1};
        end
    else
       nucs=1;
    end

    geff_sel=EPR("geff_sel");
    B_sel=EPR("B_sel");
    HF_zz_sel=EPR("HF_zz_sel");
    HF_zy_sel=EPR("HF_zy_sel");
    HF_zx_sel=EPR("HF_zx_sel");

    NQI_zz_sel=EPR("NQI_zz_sel");
    NQI_sel=EPR("NQI_sel");
    CS_zz_sel=EPR("CS_zz_sel");
    D_zz_sel=EPR("D_zz_sel");

    S_sel=EPR("S_sel");


    offsets_sel=EPR("offsets");
    Npts_EN=paramsENDOR("Npts_EN");

    v_L=paramsENDOR("v_L");


    endor_amp=zeros(1,Npts_EN);
    if length(B_sel)==0
        % No resonance orientations were found
        return
    end
    warning('This is a beta version and has not been fully tested.');


    % loop to repeat the calculation for every orientation
    parfor j=1:length(B_sel)
        endor_amp_tmp=zeros(1,Npts_EN);

        % set parameters for this orientation
        geff=geff_sel(j);
        B=B_sel(j);

        HF_zz=HF_zz_sel(j,:);

        HF_zy=HF_zy_sel(j,:);
        HF_zx=HF_zx_sel(j,:);

        NQI_zz=NQI_zz_sel(j,:);
        NQI=zeros(n_endor,3,3);

        NQI(:,:,:)=2*pi*NQI_sel(j,:,:,:);
        CS_zz=CS_zz_sel(j,:)*1e-12;
        D_zz=D_zz_sel(j,:);

        S=S_sel(j);
        offsets=offsets_sel(j,:);


        % loop over all offsets (spin manifolds)
        for i=1 : size(offsets,2)

            v_off_S=offsets(i);
            off_1=offsets(1);

            rho0=kehl_rho0(constants,paramsENDOR,B,geff,operator_spin_system,parameters,HF_zz,HF_zy,HF_zx,NQI_zz);

            %electron T1
            RT1e=kehl_relax_t1(rho0,Sx_D,parameters.T1e);
            RT1n=zeros(size(RT1e));
            RT2e=kehl_relax_t2(Sx_D,parameters.T2e);
            RT2n=zeros(size(RT2e));

            if n_spin_systems==1
                for mn=1:n_endor

                   %nuclear T1
                   RT1n=RT1e+kehl_relax_t1(rho0,Ix_D{mn},parameters.T1n);
                   RT2e=RT2e+kehl_relax_t2(Sx_D*Ix_D{mn},parameters.T2dq);
                   RT2n=RT2n+kehl_relax_t2(Ix_D{mn},parameters.T2n);
                end
            else
                RT1n=RT1e+kehl_relax_t1(rho0,Ix_D{1},parameters.T1n);
                RT2e=RT2e+kehl_relax_t2(Sx_D*Ix_D{1},parameters.T2dq);
                RT2n=RT2n+kehl_relax_t2(Ix_D{1},parameters.T2n);
            end
            R=RT1e+RT1n+RT2e+RT2n;


            start_EN=paramsENDOR("start_EN");
            step_EN=paramsENDOR("step_EN");

            oneE=expt("oneE");
            oneN=expt("oneN");

            Hfree_p=2*pi*v_off_S*Sz;
            if n_spin_systems>1
                nuc=round(i/(2*I(1)+1)) ;
                % HF
                Hfree_p=Hfree_p-2*pi*v_L(1)*Iz{1}+2*pi*v_L(1)*Iz{1}*CS_zz(nuc)+2*pi*HF_zz(nuc)*(Sz*Iz{1})+2*pi*HF_zy(nuc)*(Sz*Iy{1})+2*pi*HF_zx(nuc)*(Sz*Ix{1});
                % NQI
                if parameters.Bterm==true
                    Hfree_p=Hfree_p+NQI(1,1,1)*Ix{1}*Ix{1};
                    Hfree_p=Hfree_p+NQI(1,1,2)*Ix{1}*Iy{1};
                    Hfree_p=Hfree_p+NQI(1,1,3)*Ix{1}*Iz{1};
                    Hfree_p=Hfree_p+NQI(1,2,1)*Iy{1}*Ix{1};
                    Hfree_p=Hfree_p+NQI(1,2,2)*Iy{1}*Iy{1};
                    Hfree_p=Hfree_p+NQI(1,2,3)*Iy{1}*Iz{1};
                    Hfree_p=Hfree_p+NQI(1,3,1)*Iz{1}*Ix{1};
                    Hfree_p=Hfree_p+NQI(1,3,2)*Iz{1}*Iy{1};
                    Hfree_p=Hfree_p+NQI(1,3,3)*Iz{1}*Iz{1};
                else
                    Hfree_p=Hfree_p+pi*NQI_zz(1)*(3*Iz{1}*Iz{1}-I(1)*(I(1)+1)*eye(size(Hfree_p)));
                end
            else
                for mn=1:n_endor
                    % HF
                    Hfree_p=Hfree_p-2*pi*v_L(mn)*Iz{mn}+2*pi*HF_zz(mn)*(Sz*Iz{mn})+2*pi*HF_zy(mn)*(Sz*Iy{mn})+2*pi*HF_zx(mn)*(Sz*Ix{mn})+2*pi*v_L(mn)*CS_zz(mn)*Iz{mn};
                    % NQI
                    if parameters.Bterm==true
                        Hfree_p=Hfree_p+NQI(mn,1,1)*Ix{mn}*Ix{mn};
                        Hfree_p=Hfree_p+NQI(mn,1,2)*Ix{mn}*Iy{mn};
                        Hfree_p=Hfree_p+NQI(mn,1,3)*Ix{mn}*Iz{mn};
                        Hfree_p=Hfree_p+NQI(mn,2,1)*Iy{mn}*Ix{mn};
                        Hfree_p=Hfree_p+NQI(mn,2,2)*Iy{mn}*Iy{mn};
                        Hfree_p=Hfree_p+NQI(mn,2,3)*Iy{mn}*Iz{mn};
                        Hfree_p=Hfree_p+NQI(mn,3,1)*Iz{mn}*Ix{mn};
                        Hfree_p=Hfree_p+NQI(mn,3,2)*Iz{mn}*Iy{mn};
                        Hfree_p=Hfree_p+NQI(mn,3,3)*Iz{mn}*Iz{mn};
                    else
                        Hfree_p=Hfree_p+pi*NQI_zz(mn)*(3*Iz{mn}*Iz{mn}-I(mn)*(I(mn)+1)*eye(size(Hfree_p)));
                    end
                end

                if parameters.dipolar_active==true
                    for mn=2:(size(D_zz,2)+1)
                        dipC=2*pi*D_zz(mn-1);

                        HD=zeros(size(Hfree_p));

                        HD=HD+Ix{1}*Ix{mn};
                        HD=HD+Ix{1}*Iy{mn};
                        HD=HD+Ix{1}*Iz{mn};
                        HD=HD+Iy{1}*Ix{mn};
                        HD=HD+Iy{1}*Iy{mn};
                        HD=HD+Iy{1}*Iz{mn};
                        HD=HD+Iz{1}*Ix{mn};
                        HD=HD+Iz{1}*Iy{mn};
                        HD=HD+Iz{1}*Iz{mn};

                        HDip=dipC*(3*Iz{1}*Iz{mn}-HD)/2;
                        Hfree_p=Hfree_p+HDip;
                    end
                end
            end

            % Integration step for the Signal to account for oscillation
            if v_off_S==0
                t11=1/(off_1*Nint);
            else
                t11=1/(v_off_S*Nint);
            end

            % Radiofrequency
            v_RF=v_L(1);

                Hcorr=zeros(size(Hfree_p));
                HRF=Hfree_p;

                if parameters.Bterm==false
                    if n_spin_systems>1
                        m=1;
                        Hcorr=Hcorr+2*pi*v_RF*Iz{m};

                        HRF=HRF+2*pi*v_RF*Iz{m}+oneN*Iy{m};
                    else
                        for mn=1:n_endor
                            Hcorr=Hcorr+2*pi*v_RF*Iz{mn};

                            HRF=HRF+2*pi*v_RF*Iz{mn}+oneN*Iy{mn};
                        end
                    end
                end

                Hfree=Hfree_p+Hcorr;

                Hnonsel=full(hilb2liouv(sparse(Hfree+oneE*Sx),'comm'));

                Hfree=full(hilb2liouv(sparse(Hfree),'comm'));



                if parameters.Bterm==false

                    U5=full(propagator(operator_spin_system,1i*sparse(R-1i*full(hilb2liouv(sparse(HRF),'comm'))),t(5)));
                    if t(5)==t(7)
                        U7=U5;
                    else
                        U7=full(propagator(operator_spin_system,1i*sparse(R-1i*full(hilb2liouv(sparse(HRF),'comm'))),t(7)));
                    end

                    U1=full(propagator(operator_spin_system,1i*sparse(R-1i*Hnonsel),t(1)));

                    if t(1)==t(3)
                        U3=U1;
                    else
                        U3=full(propagator(operator_spin_system,1i*sparse(R-1i*Hnonsel),t(3)));
                    end

                    if t(1)==t(9)
                        U9=U1;
                    elseif t(3)==t(9)
                        U9=U3;
                    else
                        U9=full(propagator(operator_spin_system,1i*sparse(R-1i*Hnonsel),t(9)));
                    end

                    U2=full(propagator(operator_spin_system,1i*sparse(R-1i*Hfree),t(2)));
                    U4=full(propagator(operator_spin_system,1i*sparse(R-1i*Hfree),t(4)));

                    if t(8)==t(4)
                        U8=U4;
                    else
                        U8=full(propagator(operator_spin_system,1i*sparse(R-1i*Hfree),t(8)));
                    end
                    U10=full(propagator(operator_spin_system,1i*sparse(R-1i*Hfree),t(2)+t(3)/2));
                else
                    U5=kehl_rf_bterm_rlx(parameters,expt,v_RF,Hfree_p,Iy,t(5),n_endor,n_spin_systems,R,operator_spin_system);
                    if t(5)==t(7)
                        U7=U5;
                    else
                        U7=kehl_rf_bterm_rlx(parameters,expt,v_RF,Hfree_p,Iy,t(7),n_endor,n_spin_systems,R,operator_spin_system);
                    end

                    U1=full(propagator(operator_spin_system,1i*sparse(R-1i*Hnonsel),t(1)));

                    if t(1)==t(3)
                        U3=U1;
                    else
                        U3=full(propagator(operator_spin_system,1i*sparse(R-1i*Hnonsel),t(3)));

                    end

                    if t(1)==t(9)
                        U9=U1;
                    elseif t(3)==t(9)
                        U9=U3;
                    else
                        U9=full(propagator(operator_spin_system,1i*sparse(R-1i*Hnonsel),t(9)));
                    end

                    U2=full(propagator(operator_spin_system,1i*sparse(R-1i*Hfree),t(2)));
                    U4=full(propagator(operator_spin_system,1i*sparse(R-1i*Hfree),t(4)));
                    if t(8)==t(4)
                        U8=U4;
                    else
                        U8=full(propagator(operator_spin_system,1i*sparse(R-1i*Hfree),t(8)));
                    end

                    U10=full(propagator(operator_spin_system,1i*sparse(R-1i*Hfree),t(2)+t(3)/2));
                end

                % Evolve the densitymatrix
                rho=hilb2liouv(rho0,'statevec');

                rho=U1*rho;
                rho=U2*rho;
                rho=U3*rho;
                rho=U4*rho;
                rho=U5*rho;
                rho_t=rho;


            % loop over different separation of the RF pulses (x-axis)
            for a=1:Npts_EN
                t6=t(6)+step_EN*(a-1);
                U6=full(propagator(operator_spin_system,1i*sparse(R-1i*Hfree),t6));

                rho=rho_t;
                rho=U6*rho;
                rho=U7*rho;
                rho=U8*rho;
                rho=U9*rho;
                rho=U9*rho;
                rho=U10*rho;


                value_Sy=0;
                for b=1
                   rho_f=reshape(rho,sqrt(size(rho,1)),sqrt(size(rho,1)));
                   value_Sy=value_Sy+(real(trace(rho_f*Sy)));
                end

                if abs(HF_zz(1))>1/parameters.T2e*0.1
                    endor_amp_tmp(a)=endor_amp_tmp(a)+(value_Sy*S/(Nint*size(offsets,2)));
                end

            end
        end
        endor_amp=endor_amp+endor_amp_tmp;
    end

end

function kehl_time_rlx_grumble(spin_system,parameters)
if (~isstruct(spin_system))||(~isfield(spin_system,'bas'))||(~isfield(spin_system,'comp'))
    error('spin_system must be a Spinach spin system structure.');
end
if ~isstruct(parameters)
    error('parameters must be a structure.');
end
end

function endor_amp=kehl_time_calc(spin_system,parameters)

    % Check consistency
    kehl_time_calc_grumble(spin_system,parameters);

    % Unpack context data
    constants=parameters.constants;
    expt=parameters.expt;
    paramsENDOR=parameters.paramsENDOR;
    EPR=parameters.epr;
    operator_spin_system=parameters.operator_spin_system;
    n_endor=parameters.n_endor;
    n_spin_systems=parameters.n_spin_systems;
    I=parameters.endor_spin_numbers;

    % Obtain electron operators directly from Spinach
    Sx=full(operator(operator_spin_system,'Lx',1));
    Sy=full(operator(operator_spin_system,'Ly',1));
    Sz=full(operator(operator_spin_system,'Lz',1));

    % Obtain nuclear operators directly from Spinach
    Ix=cell(1,operator_spin_system.comp.nspins-1);
    Iy=cell(1,operator_spin_system.comp.nspins-1);
    Iz=cell(1,operator_spin_system.comp.nspins-1);
    for n=1:numel(Ix)
        spin_idx=n+1;
        Ix{n}=full(operator(operator_spin_system,'Lx',spin_idx));
        Iy{n}=full(operator(operator_spin_system,'Ly',spin_idx));
        Iz{n}=full(operator(operator_spin_system,'Lz',spin_idx));
    end
    t=expt("t");
    Nint=8;

    if n_spin_systems>1
       nucs=n_endor;
        for i=1:n_endor-1
            Ix{i+1}=Ix{1};
            Iy{i+1}=Iy{1};
            Iz{i+1}=Iz{1};
        end
    else
       nucs=1;
    end


    geff_sel=EPR("geff_sel");
    B_sel=EPR("B_sel");
    HF_zz_sel=EPR("HF_zz_sel");
    HF_zy_sel=EPR("HF_zy_sel");
    HF_zx_sel=EPR("HF_zx_sel");

    NQI_zz_sel=EPR("NQI_zz_sel");
    NQI_sel=EPR("NQI_sel");
    CS_zz_sel=EPR("CS_zz_sel");
    D_zz_sel=EPR("D_zz_sel");

    if parameters.cs_active
       [v_cs,d_cs]=eig(parameters.cs_matrix);
       cs_iso=trace(d_cs)/3*1e-12;
    else
        cs_iso=0;
    end


    S_sel=EPR("S_sel");


    offsets_sel=EPR("offsets");
    Npts_EN=paramsENDOR("Npts_EN");

    v_L=paramsENDOR("v_L");


    endor_amp=zeros(1,Npts_EN);

    if length(B_sel)==0
        % No resonance orientations were found
        return
    end


    % loop to repeat the calculation for every orientation
    parfor j=1:length(B_sel)
        % for j=1
        % set parameters for this orientation
        endor_amp_tmp=zeros(1,Npts_EN);

        geff=geff_sel(j);
        B=B_sel(j);

        HF_zz=HF_zz_sel(j,:);
        HF_zy=HF_zy_sel(j,:);
        HF_zx=HF_zx_sel(j,:);

        NQI_zz=NQI_zz_sel(j,:);
        NQI=zeros(n_endor,3,3);

        NQI(:,:,:)=2*pi*NQI_sel(j,:,:,:);

        CS_zz=CS_zz_sel(j,:)*1e-12;
        D_zz=D_zz_sel(j,:);


        S=S_sel(j);
        offsets=offsets_sel(j,:);


        % loop over all offsets (spin manifolds)
        for i=1 : size(offsets,2)

            v_off_S=offsets(i);
            off_1=offsets(1);

            [rho0]=kehl_rho0(constants,paramsENDOR,B,geff,operator_spin_system,parameters,HF_zz,HF_zy,HF_zx,NQI_zz);

            start_EN=paramsENDOR("start_EN");
            step_EN=paramsENDOR("step_EN");

            oneE=expt("oneE");
            oneN=expt("oneN");

            Hfree_p=2*pi*v_off_S*Sz;
            if n_spin_systems>1
                nuc=round(i/(2*I(1)+1)) ;
                % HF
                Hfree_p=Hfree_p-2*pi*v_L(1)*Iz{1}+2*pi*v_L(1)*Iz{1}*CS_zz(nuc)+2*pi*HF_zz(nuc)*(Sz*Iz{1})+2*pi*HF_zy(nuc)*(Sz*Iy{1})+2*pi*HF_zx(nuc)*(Sz*Ix{1});
                % NQI
                if parameters.Bterm==true
                    Hfree_p=Hfree_p+NQI(1,1,1)*Ix{1}*Ix{1};
                    Hfree_p=Hfree_p+NQI(1,1,2)*Ix{1}*Iy{1};
                    Hfree_p=Hfree_p+NQI(1,1,3)*Ix{1}*Iz{1};
                    Hfree_p=Hfree_p+NQI(1,2,1)*Iy{1}*Ix{1};
                    Hfree_p=Hfree_p+NQI(1,2,2)*Iy{1}*Iy{1};
                    Hfree_p=Hfree_p+NQI(1,2,3)*Iy{1}*Iz{1};
                    Hfree_p=Hfree_p+NQI(1,3,1)*Iz{1}*Ix{1};
                    Hfree_p=Hfree_p+NQI(1,3,2)*Iz{1}*Iy{1};
                    Hfree_p=Hfree_p+NQI(1,3,3)*Iz{1}*Iz{1};
                else
                    Hfree_p=Hfree_p+pi*NQI_zz(1)*(3*Iz{1}*Iz{1}-I(1)*(I(1)+1)*eye(size(Hfree_p)));
                end
            else
                for mm=1:n_endor
                    % HF
                    Hfree_p=Hfree_p-2*pi*v_L(mm)*Iz{mm}+2*pi*HF_zz(mm)*(Sz*Iz{mm})+2*pi*HF_zy(mm)*(Sz*Iy{mm})+2*pi*HF_zx(mm)*(Sz*Ix{mm})+2*pi*v_L(mm)*CS_zz(mm)*Iz{mm};
                    % NQI
                    if parameters.Bterm==true
                        HNQI=zeros(size(Hfree_p));
                        HNQI=HNQI+Ix{mm}*NQI(mm,1,1)*Ix{mm};
                        HNQI=HNQI+Ix{mm}*NQI(mm,1,2)*Iy{mm};
                        HNQI=HNQI+Ix{mm}*NQI(mm,1,3)*Iz{mm};
                        HNQI=HNQI+Iy{mm}*NQI(mm,2,1)*Ix{mm};
                        HNQI=HNQI+Iy{mm}*NQI(mm,2,2)*Iy{mm};
                        HNQI=HNQI+Iy{mm}*NQI(mm,2,3)*Iz{mm};
                        HNQI=HNQI+Iz{mm}*NQI(mm,3,1)*Ix{mm};
                        HNQI=HNQI+Iz{mm}*NQI(mm,3,2)*Iy{mm};
                        HNQI=HNQI+Iz{mm}*NQI(mm,3,3)*Iz{mm};
                        Hfree_p=Hfree_p+HNQI;
                    else
                        Hfree_p=Hfree_p+pi*NQI_zz(mm)*(3*Iz{mm}*Iz{mm}-I(mm)*(I(mm)+1)*eye(size(Hfree_p)));

                    end

                end

                if parameters.dipolar_active==true
                    for mm=2:(size(D_zz,2)+1)
                        dipC=2*pi*D_zz(mm-1);

                        HD=zeros(size(Hfree_p));

                        HD=HD+Ix{1}*Ix{mm};
                        HD=HD+Ix{1}*Iy{mm};
                        HD=HD+Ix{1}*Iz{mm};
                        HD=HD+Iy{1}*Ix{mm};
                        HD=HD+Iy{1}*Iy{mm};
                        HD=HD+Iy{1}*Iz{mm};
                        HD=HD+Iz{1}*Ix{mm};
                        HD=HD+Iz{1}*Iy{mm};
                        HD=HD+Iz{1}*Iz{mm};

                        HDip=dipC*(3*Iz{1}*Iz{mm}-HD)/2;

                        Hfree_p=Hfree_p+HDip;
                    end
                end
            end
                v_RF=v_L(1);
                Hcorr=zeros(size(Hfree_p));
                if n_spin_systems>1
                    m=1;
                    Hcorr=Hcorr+2*pi*v_RF*Iz{m};
                else

                    for mm=1:n_endor
                        Hcorr=Hcorr+2*pi*v_RF*Iz{mm};

                    end
                end

                if parameters.Bterm==false
                    Hfree=Hfree_p+Hcorr;
                else
                    Hfree=Hfree_p;
                end

                % mw pulses
                Hnonsel=Hfree+oneE*Sx;

                HRF=Hfree;

                if parameters.Bterm==false
                    if n_spin_systems>1
                        m=1;
                        HRF=HRF+oneN*Iy{m};
                    else
                        for mm=1:n_endor
                            HRF=HRF+oneN*Iy{mm};
                        end
                    end
                end


                % Integration step for the Signal to account for oscillation
                if v_off_S==0
                    t11=1/abs(off_1*Nint);
                else
                    t11=1/abs(v_off_S*Nint);
                end


                if parameters.Bterm==false
                    U5=full(propagator(operator_spin_system,sparse(HRF),t(5)));
                    if t(5)==t(7)
                        U7=U5;
                    else
                        U7=full(propagator(operator_spin_system,sparse(HRF),t(7)));
                    end

                else
                    U5=kehl_rf_bterm(parameters,expt,v_RF,Hfree,Iy,t(5),n_endor,n_spin_systems,operator_spin_system);
                    if t(5)==t(7)
                        U7=U5;
                    else
                        U7=kehl_rf_bterm(parameters,expt,v_RF,Hfree,Iy,t(7),n_endor,n_spin_systems,operator_spin_system);
                    end
                end

                % Calculate the propagators
                U1=full(propagator(operator_spin_system,sparse(Hnonsel),t(1)));

                if t(1)==t(3)
                    U3=U1;
                else
                    U3=full(propagator(operator_spin_system,sparse(Hnonsel),t(3)));
                end

                if t(1)==t(9)
                    U9=U1;
                elseif t(3)==t(9)
                    U9=U3;
                else
                    U9=full(propagator(operator_spin_system,sparse(Hfree),t(9)));
                end

                U2=full(propagator(operator_spin_system,sparse(Hfree),t(2)));

                U4=full(propagator(operator_spin_system,sparse(Hfree),t(4)));

                U10=full(propagator(operator_spin_system,sparse(Hfree),t(2)+t(3)/2));

                U11=full(propagator(operator_spin_system,sparse(Hfree),t11));


                % Evolve the densitymatrix
                rho=rho0;
                rho=U1*rho*U1';
                rho=U2*rho*U2';
                rho=U3*rho*U3';
                rho=U4*rho*U4';
                rho=U5*rho*U5';

                if parameters.Bterm==true
                    rho=diag(diag(rho));
                end

            rho_t=rho;

             % loop over different separation of the RF pulses (x-axis)
            for a=1:Npts_EN
                t6=t(6)+step_EN*(a-1);
                t8=t(8)-step_EN*(a-1);


                if parameters.Bterm==false
                    U6=full(propagator(operator_spin_system,sparse(Hfree_p+Hcorr),t6));
                else
                    Hfree=Hfree_p;
                    U6=full(propagator(operator_spin_system,sparse(Hfree),t6));
                end
                U8=full(propagator(operator_spin_system,sparse(Hfree),t8));


                rho=rho_t;
                rho=U6*rho*U6';
                rho=U7*rho*U7';

                if parameters.Bterm==true
                    rho=diag(diag(rho));
                end

                rho=U8*rho*U8';
                rho=U9*rho*U9';
                rho=U10*rho*U10';

                value_Sy=0;
                for b=1:Nint
                   rho=U11*rho*U11';
                   value_Sy=value_Sy+(real(trace(rho*Sy)));
                end

                endor_amp_tmp(a)=endor_amp_tmp(a)+(value_Sy*S/(Nint*size(offsets,2)));
            end
        end
        endor_amp=endor_amp+endor_amp_tmp;
    end

end

function kehl_time_calc_grumble(spin_system,parameters)
if (~isstruct(spin_system))||(~isfield(spin_system,'bas'))||(~isfield(spin_system,'comp'))
    error('spin_system must be a Spinach spin system structure.');
end
if ~isstruct(parameters)
    error('parameters must be a structure.');
end
end
