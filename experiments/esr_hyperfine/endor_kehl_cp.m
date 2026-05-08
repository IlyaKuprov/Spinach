%ENDOR_KEHL_CP cross-polarisation ENDOR pulse sequence for the Kehl ENDOR context.
%
%   ENDOR_AMP=ENDOR_KEHL_CP(SPIN_SYSTEM,PARAMETERS,H,R,K) is called by
%   endor_kehl_context.m. H, R, and K are present for Spinach experiment
%   signature compatibility; this sequence uses orientation-selected
%   effective Hamiltonian data prepared by the context.


function endor_amp=endor_kehl_cp(spin_system,parameters,H,R,K)

% Check consistency
grumble(spin_system,parameters,H,R,K);
if parameters.Relax==true
    endor_amp=kehl_cp_calc_rlx(spin_system,parameters);
else
    endor_amp=kehl_cp_calc(spin_system,parameters);
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

function endor_amp=kehl_cp_calc_rlx(spin_system,parameters)

    % Check consistency
    kehl_cp_calc_rlx_grumble(spin_system,parameters);

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

    % get parameters from Maps



    t=expt("t");
    Nint=8;

    n_spin_systems=n_spin_systems;

    n_endor=n_endor;
    Npts_CP=expt("Npts_CP");

    if expt("Npts_CP")>1
        step_CP=expt("range_CP")/(expt("Npts_CP")-1);
    else
        step_CP=0;
    end

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
    % length(B_sel)

    if length(B_sel)==0
        % No resonance orientations were found
        return
    end

    % loop to repeat the calculation for every orientation
    for j=1:length(B_sel)
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

            % % loop over nuclei if not all in one spin system
        % for n = 1 : nucs
        % loop over all offsets (spin manifolds)
        for i=1 : size(offsets,2)

            v_off_S=offsets(i);
            off_1=offsets(1);


            rho0=kehl_rho0(constants,paramsENDOR,B,geff,operator_spin_system,parameters,HF_zz,HF_zy,HF_zx,NQI_zz);

            %electron T1
            RT1e=kehl_relax_t1(rho0,Sx_D,parameters.T1e);
            RT1n=zeros(size(RT1e));

            % electron T2
            RT2e=kehl_relax_t2(Sx_D,parameters.T2e);
            RT2n=zeros(size(RT2e));

            if n_spin_systems==1
                for mm=1:n_endor

                   %nuclear T1
                   RT1n=RT1e+kehl_relax_t1(rho0,Ix_D{mm},parameters.T1n);

                   % double and zero quantum T2
                   RT2e=RT2e+kehl_relax_t2(Sx_D*Ix_D{mm},parameters.T2dq);

                   % nuclear T2
                   RT2n=RT2n+kehl_relax_t2(Ix_D{mm},parameters.T2n);
                end
            else

                %nuclear T1
                RT1n=RT1e+kehl_relax_t1(rho0,Ix_D{1},parameters.T1n);

                % double and zero quantum T2
                RT2e=RT2e+kehl_relax_t2(Sx_D*Ix_D{1},parameters.T2dq);

                % nuclear T2
                RT2n=RT2n+kehl_relax_t2(Ix_D{1},parameters.T2n);
            end

            % full relaxation superoperator
            R=RT1e+RT1n+RT2e+RT2n;

            % for tilted frame Rho relaxation times

            %electron T1
            RT1eRho=kehl_relax_t1(rho0,Sx_D,parameters.T1eR);
            RT1nRho=zeros(size(RT1eRho));

            % electron T2
            RT2eRho=kehl_relax_t2(Sx_D,parameters.T2eR);
            RT2nRho=zeros(size(RT2eRho));

            if n_spin_systems==1
                for mm=1:n_endor

                   %nuclear T1
                   RT1nRho=RT1eRho+kehl_relax_t1(rho0,Ix_D{mm},parameters.T1nR);

                   % double and zero quantum T2
                   RT2eRho=RT2eRho+kehl_relax_t2(Sx_D*Ix_D{mm},parameters.T2dqR);

                   % nuclear T2
                   RT2nRho=RT2nRho+kehl_relax_t2(Ix_D{mm},parameters.T2nR);
                end
            else

                %nuclear T1
                RT1nRho=RT1eRho+kehl_relax_t1(rho0,Ix_D{1},parameters.T1nR);

                % double and zero quantum T2
                RT2eRho=RT2eRho+kehl_relax_t2(Sx_D*Ix_D{1},parameters.T2dqR);

                % nuclear T2
                RT2nRho=RT2nRho+kehl_relax_t2(Ix_D{1},parameters.T2nR);
            end

            % full relaxation superoperator
            RRho=RT1eRho+RT1nRho+RT2eRho+RT2nRho;

            start_EN=paramsENDOR("start_EN");
            step_EN=paramsENDOR("step_EN");

            prep=expt("prep");
            sl=expt("SL");
            cp=expt("CP");
            oneE=expt("oneE");
            oneN=expt("oneN");

            Hfree_p=2*pi*v_off_S*Sz;
            if n_spin_systems>1
                if parameters.powder==1
                    mn=round(i/(2*I(1)+1));
                else
                    mn=i;
                end
                % HF
                Hfree_p=Hfree_p-2*pi*v_L(mn)*Iz{1}+2*pi*v_L(mn)*Iz{1}*CS_zz(mn)+2*pi*HF_zz(mn)*(Sz*Iz{1})+2*pi*HF_zy(mn)*(Sz*Iy{1})+2*pi*HF_zx(mn)*(Sz*Ix{1});
                % NQI
                if parameters.Bterm==true
                    Hfree_p=Hfree_p+NQI(mn,1,1)*Ix{1}*Ix{1};
                    Hfree_p=Hfree_p+NQI(mn,1,2)*Ix{1}*Iy{1};
                    Hfree_p=Hfree_p+NQI(mn,1,3)*Ix{1}*Iz{1};
                    Hfree_p=Hfree_p+NQI(mn,2,1)*Iy{1}*Ix{1};
                    Hfree_p=Hfree_p+NQI(mn,2,2)*Iy{1}*Iy{1};
                    Hfree_p=Hfree_p+NQI(mn,2,3)*Iy{1}*Iz{1};
                    Hfree_p=Hfree_p+NQI(mn,3,1)*Iz{1}*Ix{1};
                    Hfree_p=Hfree_p+NQI(mn,3,2)*Iz{1}*Iy{1};
                    Hfree_p=Hfree_p+NQI(mn,3,3)*Iz{1}*Iz{1};
                else
                    Hfree_p=Hfree_p+pi*NQI_zz(mn)*(3*Iz{1}*Iz{1}-I(1)*(I(1)+1)*eye(size(Hfree_p)));
                end

                % get offset for sc CP calculation
                if parameters.powder==false
                    [v_CP,paramsENDOR]=kehl_cp_offset(parameters,paramsENDOR,expt,v_off_S,HF_zz,NQI_zz,mn,i);
                end

            else
                for mm=1:n_endor
                    % HF
                    Hfree_p=Hfree_p-2*pi*v_L(mm)*Iz{mm}+2*pi*HF_zz(mm)*(Sz*Iz{mm})+2*pi*HF_zy(mm)*(Sz*Iy{mm})+2*pi*HF_zx(mm)*(Sz*Ix{mm})+2*pi*v_L(mm)*CS_zz(mm)*Iz{mm};
                    % NQI
                    if parameters.Bterm==true
                        Hfree_p=Hfree_p+NQI(mm,1,1)*Ix{mm}*Ix{mm};
                        Hfree_p=Hfree_p+NQI(mm,1,2)*Ix{mm}*Iy{mm};
                        Hfree_p=Hfree_p+NQI(mm,1,3)*Ix{mm}*Iz{mm};
                        Hfree_p=Hfree_p+NQI(mm,2,1)*Iy{mm}*Ix{mm};
                        Hfree_p=Hfree_p+NQI(mm,2,2)*Iy{mm}*Iy{mm};
                        Hfree_p=Hfree_p+NQI(mm,2,3)*Iy{mm}*Iz{mm};
                        Hfree_p=Hfree_p+NQI(mm,3,1)*Iz{mm}*Ix{mm};
                        Hfree_p=Hfree_p+NQI(mm,3,2)*Iz{mm}*Iy{mm};
                        Hfree_p=Hfree_p+NQI(mm,3,3)*Iz{mm}*Iz{mm};
                    else
                        Hfree_p=Hfree_p+pi*NQI_zz(mm)*(3*Iz{mm}*Iz{mm}-I(mm)*(I(mm)+1)*eye(size(Hfree_p)));
                    end

                    % get offset for sc CP calculation
                    if parameters.powder==false
                        [v_CP,paramsENDOR]=kehl_cp_offset(parameters,paramsENDOR,expt,v_off_S,HF_zz,NQI_zz,mm,i);
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

            % Integration step for the Signal to account for oscillation
            if v_off_S==0
                t11=1/(off_1*Nint);
            else
                t11=1/(v_off_S*Nint);
            end

            % loop over cp frequencies (y-axis)
            for c=1:Npts_CP

            if parameters.powder==true
                v_CP=expt("start_CP")+step_CP*(c-1);
            end

            % loop over rf frequencies (x-axis)
            parfor a=1:Npts_EN
            % for a=1   %%% for testing

                if abs(HF_zz(1))>1/parameters.T2e*0.1

                % Radiofrequency
                v_RF=(start_EN+step_EN*(a-1));

                Hcorr=zeros(size(Hfree_p));
                % Hamiltonian for RF pulse (no HF enhancement)
                HRF=Hfree_p;
                HRFb=zeros(size(Hfree_p));
                HSL=Hfree_p+sl*Sy;


                if parameters.Bterm==false
                    if n_spin_systems>1
                        m=1;
                        HRF=HRF+2*pi*v_RF*Iz{m}+oneN*Ix{m};
                        HRFb=2*pi*(v_RF-v_CP)*Iz{m};
                        HSL=HSL+2*pi*v_CP*Iz{m}+cp*Ix{m};
                    else
                        for mm=1:n_endor
                            HRF=HRF+2*pi*v_RF*Iz{mm}+oneN*Ix{mm};
                            HRFb=2*pi*(v_RF-v_CP)*Iz{mm};
                            HSL=HSL+2*pi*v_CP*Iz{mm}+cp*Ix{mm};
                        end
                    end
                end

                [V,D]=eig(HSL);
                if ~issorted(diag(D))
                  [d,inds]=sort(diag(D));
                  D=diag(d);
                  V=V(:, inds);
                end
                %eigenvalues are order from negative large to positive large!
                W=V^(-1);

                HSL_t=W*HSL*W^(-1);


                if n_spin_systems>1
                    m=1;
                    Hcorr=Hcorr+2*pi*v_CP*Iz{m};
                else

                    for mm=1:n_endor
                        Hcorr=Hcorr+2*pi*v_CP*Iz{mm};
                    end
                end


                Hfree=Hfree_p+Hcorr;

                Hnonsel=full(hilb2liouv(sparse(Hfree+oneE*Sx),'comm'));
                Hprep=full(hilb2liouv(sparse(Hfree+prep*Sx),'comm'));

                Hfree=full(hilb2liouv(sparse(Hfree),'comm'));


                    if parameters.Bterm==false
                        U3=full(propagator(operator_spin_system,1i*sparse(RRho-1i*full(hilb2liouv(sparse(HSL_t),'comm'))),t(3)));

                        U5=full(propagator(operator_spin_system,1i*sparse(R-1i*full(hilb2liouv(sparse(HRF),'comm'))),t(5)));

                        U1=full(propagator(operator_spin_system,1i*sparse(R-1i*Hprep),t(1)));
                        U2=full(propagator(operator_spin_system,1i*sparse(R-1i*Hfree),t(2)));

                        U4=full(propagator(operator_spin_system,1i*sparse(R-1i*Hfree),t(4)));

                        U6=full(propagator(operator_spin_system,1i*sparse(R-1i*Hfree),t(6)));
                        U7=full(propagator(operator_spin_system,1i*sparse(R-1i*Hnonsel),t(7)));
                        U8=full(propagator(operator_spin_system,1i*sparse(R-1i*Hfree),t(8)));
                        U9=full(propagator(operator_spin_system,1i*sparse(R-1i*Hnonsel),t(9)));
                        U10=full(propagator(operator_spin_system,1i*sparse(R-1i*Hfree),t(8)+t(7)/2));
                        U11=full(propagator(operator_spin_system,1i*sparse(R-1i*Hfree),t11));
                    else
                        t_stepCP=1/(v_CP*parameters.N_stepRF);
                        H_SL=HSL_t;
                        U_SL=eye(size(HSL_t,1)^2,size(HSL_t,1)^2);
                        for ll=1:parameters.N_stepRF
                            if n_spin_systems==1
                                for mm=1:n_endor
                                    H_SL=H_SL+2*expt('CP')*Ix{mm}*cos(2*pi*v_CP*t_stepCP*(ll-1));
                                end
                            else
                                H_SL=H_SL+2*expt('CP')*Ix{1}*cos(2*pi*v_CP*t_stepCP*(ll-1));
                            end
                            G=RRho/t_stepCP-1i*full(hilb2liouv(sparse(H_SL),'comm'));
                            U_step=full(propagator(operator_spin_system,1i*sparse(G),t_stepCP));
                            U_SL=U_step*U_SL;
                        end
                        U3=(U_SL^(t(3)*v_CP));
                        U5=kehl_rf_bterm_rlx(parameters,expt,v_RF,Hfree_p,Iy,t(5),n_endor,n_spin_systems,R,operator_spin_system);

                        U1=full(propagator(operator_spin_system,1i*sparse(R-1i*Hprep),t(1)));
                        U2=full(propagator(operator_spin_system,1i*sparse(R-1i*Hfree),t(2)));

                        U4=full(propagator(operator_spin_system,1i*sparse(R-1i*Hfree),t(4)));

                        U6=full(propagator(operator_spin_system,1i*sparse(R-1i*Hfree),t(6)));
                        U7=full(propagator(operator_spin_system,1i*sparse(R-1i*Hnonsel),t(7)));
                        U8=full(propagator(operator_spin_system,1i*sparse(R-1i*Hfree),t(8)));
                        U9=full(propagator(operator_spin_system,1i*sparse(R-1i*Hnonsel),t(9)));
                        U10=full(propagator(operator_spin_system,1i*sparse(R-1i*Hfree),t(8)+t(7)/2));
                        U11=full(propagator(operator_spin_system,1i*sparse(R-1i*Hfree),t11));

                    end

                % Evolve the densitymatrix
                rho=hilb2liouv(rho0,'statevec');
                rho=U1*rho;
                rho=U2*rho;

                rho_t=W*reshape(rho,sqrt(size(rho,1)),sqrt(size(rho,1)))*W^(-1);
                rho_t=hilb2liouv(rho_t,'statevec');
                rho_t=U3*rho_t;
                rho=V*reshape(rho_t,sqrt(size(rho_t,1)),sqrt(size(rho_t,1)))*V^(-1);
                rho=hilb2liouv(rho,'statevec');

                rho=U4*rho;
                rho=U5*rho;
                rho=U6*rho;
                rho=U7*rho;
                rho=U8*rho;
                rho=U9*rho;
                rho=U10*rho;

                value_Sy=0;
                for b=1:Nint
                   rho=U11*rho;
                   rho_f=reshape(rho,sqrt(size(rho,1)),sqrt(size(rho,1)));
                   value_Sy=value_Sy+real(trace(rho_f*Sy));
                end
                    endor_amp(a)=endor_amp(a)+(value_Sy*S/(Nint*size(offsets,2)));
                end
            end
            end

        end
    end
end

function kehl_cp_calc_rlx_grumble(spin_system,parameters)
if (~isstruct(spin_system))||(~isfield(spin_system,'bas'))||(~isfield(spin_system,'comp'))
    error('spin_system must be a Spinach spin system structure.');
end
if ~isstruct(parameters)
    error('parameters must be a structure.');
end
end

function endor_amp=kehl_cp_calc(spin_system,parameters)

    % Check consistency
    kehl_cp_calc_grumble(spin_system,parameters);

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

    n_spin_systems=n_spin_systems;

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
    Npts_CP=expt("Npts_CP");

    if expt("Npts_CP")>1
        step_CP=expt("range_CP")/(expt("Npts_CP")-1);
    else
        step_CP=0;
    end

    v_L=paramsENDOR("v_L");


    endor_amp=zeros(1,Npts_EN);
     length(B_sel)
    if length(B_sel)==0
        % No resonance orientations were found
        return
    end


    % loop to repeat the calculation for every orientation
    for j=1:length(B_sel)
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

        % % loop over nuclei if not all in one spin system
        % for n = 1:nucs

        % loop over all offsets (spin manifolds)
        for i=1 : size(offsets,2)

            v_off_S=offsets(i);
            off_1=offsets(1);

            rho0=kehl_rho0(constants,paramsENDOR,B,geff,operator_spin_system,parameters,HF_zz,HF_zy,HF_zx,NQI_zz);

            start_EN=paramsENDOR("start_EN");
            step_EN=paramsENDOR("step_EN");

            prep=expt("prep");
            sl=expt("SL");
            cp=expt("CP");
            oneE=expt("oneE");
            oneN=expt("oneN");

            % define free Hamiltonian
            Hfree_p=2*pi*v_off_S*Sz;
            if n_spin_systems>1
                if parameters.powder==1
                    mn=round(i/(2*I(1)+1));
                else
                    mn=i;
                end
                % HF
                Hfree_p=Hfree_p-2*pi*v_L(mn)*Iz{1}+2*pi*v_L(mn)*Iz{1}*CS_zz(mn)+2*pi*HF_zz(mn)*(Sz*Iz{1})+2*pi*HF_zy(mn)*(Sz*Iy{1})+2*pi*HF_zx(mn)*(Sz*Ix{1});
                % NQI
                if parameters.Bterm==true
                    Hfree_p=Hfree_p+NQI(mn,1,1)*Ix{1}*Ix{1};
                    Hfree_p=Hfree_p+NQI(mn,1,2)*Ix{1}*Iy{1};
                    Hfree_p=Hfree_p+NQI(mn,1,3)*Ix{1}*Iz{1};
                    Hfree_p=Hfree_p+NQI(mn,2,1)*Iy{1}*Ix{1};
                    Hfree_p=Hfree_p+NQI(mn,2,2)*Iy{1}*Iy{1};
                    Hfree_p=Hfree_p+NQI(mn,2,3)*Iy{1}*Iz{1};
                    Hfree_p=Hfree_p+NQI(mn,3,1)*Iz{1}*Ix{1};
                    Hfree_p=Hfree_p+NQI(mn,3,2)*Iz{1}*Iy{1};
                    Hfree_p=Hfree_p+NQI(mn,3,3)*Iz{1}*Iz{1};
                else
                    Hfree_p=Hfree_p+pi*NQI_zz(mn)*(3*Iz{1}*Iz{1}-I(1)*(I(1)+1)*eye(size(Hfree_p)));
                end

                % get offset for sc CP calculation
                if parameters.powder==false
                    [v_CP,paramsENDOR]=kehl_cp_offset(parameters,paramsENDOR,expt,v_off_S,HF_zz,NQI_zz,mn,i);
                end
            else
                for mm=1:n_endor
                    % HF
                    Hfree_p=Hfree_p-2*pi*v_L(mm)*Iz{mm}+2*pi*HF_zz(mm)*(Sz*Iz{mm})+2*pi*HF_zy(mm)*(Sz*Iy{mm})+2*pi*HF_zx(mm)*(Sz*Ix{mm})+2*pi*v_L(mm)*CS_zz(mm)*Iz{mm};
                    % NQI
                    if parameters.Bterm==true
                        Hfree_p=Hfree_p+NQI(mm,1,1)*Ix{mm}*Ix{mm};
                        Hfree_p=Hfree_p+NQI(mm,1,2)*Ix{mm}*Iy{mm};
                        Hfree_p=Hfree_p+NQI(mm,1,3)*Ix{mm}*Iz{mm};
                        Hfree_p=Hfree_p+NQI(mm,2,1)*Iy{mm}*Ix{mm};
                        Hfree_p=Hfree_p+NQI(mm,2,2)*Iy{mm}*Iy{mm};
                        Hfree_p=Hfree_p+NQI(mm,2,3)*Iy{mm}*Iz{mm};
                        Hfree_p=Hfree_p+NQI(mm,3,1)*Iz{mm}*Ix{mm};
                        Hfree_p=Hfree_p+NQI(mm,3,2)*Iz{mm}*Iy{mm};
                        Hfree_p=Hfree_p+NQI(mm,3,3)*Iz{mm}*Iz{mm};
                    else
                        Hfree_p=Hfree_p+pi*NQI_zz(mm)*(3*Iz{mm}*Iz{mm}-I(mm)*(I(mm)+1)*eye(size(Hfree_p)));
                    end

                    % get offset for sc CP calculation
                    if parameters.powder==false
                        [v_CP,paramsENDOR]=kehl_cp_offset(parameters,paramsENDOR,expt,v_off_S,HF_zz,NQI_zz,mm,i);
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

                % Hamiltonian for mw pulses
                Hprep_p=Hfree_p+prep*Sx;
                Hnonsel_p=Hfree_p+oneE*Sx;

                % Integration step for the Signal to account for oscillation
                if v_off_S==0
                    t11=1/(off_1*Nint);
                else
                    t11=1/(v_off_S*Nint);
                end


                % Calculate the propagators
                U1_p=full(propagator(operator_spin_system,sparse(Hprep_p),t(1)));
                U2_p=full(propagator(operator_spin_system,sparse(Hfree_p),t(2)));
                % SL/CP step
                U4_p=full(propagator(operator_spin_system,sparse(Hfree_p),t(4)));
                % RF pulse
                U6_p=full(propagator(operator_spin_system,sparse(Hfree_p),t(6)));
                U7_p=full(propagator(operator_spin_system,sparse(Hnonsel_p),t(7)));
                U8_p=full(propagator(operator_spin_system,sparse(Hfree_p),t(8)));
                U9_p=full(propagator(operator_spin_system,sparse(Hnonsel_p),t(9)));
                U10_p=full(propagator(operator_spin_system,sparse(Hfree_p),t(8)+t(7)/2));
                U11_p=full(propagator(operator_spin_system,sparse(Hfree_p),t11));

            % loop over cp frequencies (y-axis)
            for c=1:Npts_CP

            if parameters.powder==true
                v_CP=expt("start_CP")+step_CP*(c-1);
            end


            % loop over rf frequencies (x-axis)
            parfor a=1:Npts_EN
            % for a=1   %%% for testing

                % Radiofrequency
                v_RF=(start_EN+step_EN*(a-1));

                % Calculate the Free evolution hamiltonian
                % (not including the Bterm)

                Hcorr=zeros(size(Hfree_p));
                % Hamiltonian for RF pulse (no HF enhancement)
                HRF=Hfree_p;
                HRFb=zeros(size(Hfree_p));

                HSL=Hfree_p+sl*Sy;

                [V,D]=eig(HSL);
                if ~issorted(diag(D))
                  [d,inds]=sort(diag(D));
                  D=diag(d);
                  V=V(:, inds);
                end
                %eigenvalues are order from negative large to positive large!
                W=V^(-1);


                if parameters.Bterm==false
                    if n_spin_systems>1
                        m=1;
                        HRF=HRF+2*pi*v_RF*Iz{m}+oneN*Ix{m};
                        HRFb=2*pi*(v_RF-v_CP)*Iz{m};
                        HSL=HSL+2*pi*v_CP*Iz{m}+cp*Ix{m};
                    else
                        for mm=1:n_endor
                            HRF=HRF+2*pi*v_RF*Iz{mm}+oneN*Ix{mm};
                            HRFb=2*pi*(v_RF-v_CP)*Iz{mm};
                            HSL=HSL+2*pi*v_CP*Iz{mm}+cp*Ix{mm};
                        end
                    end
                end

                if n_spin_systems>1
                    m=1;
                    Hcorr=Hcorr+2*pi*v_CP*Iz{m};
                else

                    for mm=1:n_endor
                        Hcorr=Hcorr+2*pi*v_CP*Iz{mm};
                    end
                end

                U5b=[];
                if parameters.Bterm==false
                    U3=full(propagator(operator_spin_system,sparse(HSL),t(3)));
                    U5=full(propagator(operator_spin_system,sparse(HRF),t(5)));
                    U5b=full(propagator(operator_spin_system,sparse(HRFb),t(5)));

                    U1=U1_p*full(propagator(operator_spin_system,sparse(Hcorr),t(1)));
                    U2=U2_p*full(propagator(operator_spin_system,sparse(Hcorr),t(2)));

                    U4=U4_p*full(propagator(operator_spin_system,sparse(Hcorr),t(4)));

                    U6=U6_p*full(propagator(operator_spin_system,sparse(Hcorr),t(6)));
                    U7=U7_p*full(propagator(operator_spin_system,sparse(Hcorr),t(7)));
                    U8=U8_p*full(propagator(operator_spin_system,sparse(Hcorr),t(8)));
                    U9=U9_p*full(propagator(operator_spin_system,sparse(Hcorr),t(9)));
                    U10=U10_p*full(propagator(operator_spin_system,sparse(Hcorr),t(8)+t(7)/2));
                    U11=U11_p*full(propagator(operator_spin_system,sparse(Hcorr),t11));
                else
                    Hfree=Hfree_p;
                    HSL_p=Hfree+sl*Sy;

                    t_stepCP=1/(v_RF*parameters.N_stepRF);
                    U_SL=eye(size(HSL_p));
                    H_SL=HSL_p;
                    for ll=1:parameters.N_stepRF
                        if n_spin_systems==1
                            for mm=1:n_endor
                                H_SL=H_SL+2*expt('CP')*Ix{mm}*cos(2*pi*v_RF*t_stepCP*(ll-1));
                            end
                        else
                            H_SL=H_SL+2*expt('CP')*Ix{1}*cos(2*pi*v_RF*t_stepCP*(ll-1));
                        end
                        U_step=full(propagator(operator_spin_system,sparse(H_SL),t_stepCP));
                        U_SL=U_step*U_SL;
                    end
                    U3=(U_SL^(t(3)*v_RF));
                    U5=kehl_rf_bterm(parameters,expt,v_RF,Hfree,Iy,t(5),n_endor,n_spin_systems,operator_spin_system);

                    U1=U1_p;
                    U2=U2_p;

                    U4=U4_p;

                    U6=U6_p;
                    U7=U7_p;
                    U8=U8_p;
                    U9=U9_p;
                    U10=U10_p;
                    U11=U11_p;
                end


                % Evolve the densitymatrix
                rho=rho0;
                rho=U1*rho*U1';
                rho=U2*rho*U2';
                if t(2)~=0
                    rho=diag(diag(rho));
                end

                rho=U3*rho*U3';

                rho_t=W*rho*W^(-1);
                rho_t=diag(diag(rho_t));
                rho=V*rho_t*V^(-1);

                rho=U4*rho*U4';
                rho=U5*rho*U5';

                if parameters.Bterm==true
                    rho=diag(diag(rho));
                else
                    rho=U5b*rho*U5b';
                end

                rho=U6*rho*U6';
                rho=U7*rho*U7';
                rho=U8*rho*U8';
                rho=U9*rho*U9';
                rho=U10*rho*U10';


                value_Sy=0;
                for b=1:Nint*10
                   rho=U11*rho*U11';
                   value_Sy=value_Sy+real(trace(rho*Sy));
                end

                if n_spin_systems>1
                    if parameters.powder==0
                        s=1;
                    else
                        s=size(offsets,2)/n_endor;
                    end
                else
                    s=size(offsets,2);
                end

                endor_amp(a)=endor_amp(a)+(value_Sy*S/(Nint*s));

            end
            end
        end

    end
end

function kehl_cp_calc_grumble(spin_system,parameters)
if (~isstruct(spin_system))||(~isfield(spin_system,'bas'))||(~isfield(spin_system,'comp'))
    error('spin_system must be a Spinach spin system structure.');
end
if ~isstruct(parameters)
    error('parameters must be a structure.');
end
end
