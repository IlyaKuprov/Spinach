%ENDOR_KEHL_SPINLOCK spin-lock ENDOR pulse sequence for the Kehl ENDOR context.
%
%   ENDOR_AMP=ENDOR_KEHL_SPINLOCK(SPIN_SYSTEM,PARAMETERS,H,R,K) is called by
%   endor_kehl_context.m. H, R, and K are present for Spinach experiment
%   signature compatibility; this sequence uses orientation-selected
%   effective Hamiltonian data prepared by the context.


function endor_amp=endor_kehl_spinlock(spin_system,parameters,H,R,K)

% Check consistency
grumble(spin_system,parameters,H,R,K);
if parameters.Relax==true
    endor_amp=kehl_spinlock_rlx(spin_system,parameters);
else
    endor_amp=kehl_spinlock_calc(spin_system,parameters);
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

function endor_amp=kehl_spinlock_rlx(spin_system,parameters)

    % Check consistency
    kehl_spinlock_rlx_grumble(spin_system,parameters);

    % Unpack legacy data
    constants=parameters.constants;
    spinOps=parameters.spinOps;
    spinSys=parameters.spinSys;
    expt=parameters.expt;
    paramsENDOR=parameters.paramsENDOR;
    EPR=parameters.epr;
    spinOps_D=kehl_diag_ops(spinOps,spinSys('Ni_ENDOR'));


    % get values from Maps
    Sz=spinOps("Sz");
    Sx=spinOps("Sx");
    Sy=spinOps("Sy");
    Iz=spinOps("Iz");
    Ix=spinOps("Ix");
    Iy=spinOps("Iy");


    Sz_D=spinOps_D("Sz");
    Sx_D=spinOps_D("Sx");
    Sy_D=spinOps_D("Sy");
    Iz_D=spinOps_D("Iz");
    Ix_D=spinOps_D("Ix");
    Iy_D=spinOps_D("Iy");


    t=expt("t");
    Nint=8;


    N_spinSys=spinSys("N_SpinSys");

    Ni_ENDOR=spinSys("Ni_ENDOR");

    if N_spinSys>1
       nucs=Ni_ENDOR;
        for i=1:Ni_ENDOR-1
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

    S_sel=EPR("S_sel");


    offsets_sel=EPR("offsets");
    Npts_EN=paramsENDOR("Npts_EN");

    I=spinSys("I");
    v_L=paramsENDOR("v_L");


    endor_amp=zeros(1,Npts_EN);

    if length(B_sel)==0
        % No resonance orientations were found
        return
    end

    % loop to repeat the calculation for every orientation
    parfor j=1:length(B_sel)
        % set parameters for this orientation
        geff=geff_sel(j);
        B=B_sel(j);

        HF_zz=HF_zz_sel(j,:);
        HF_zy=HF_zy_sel(j,:);
        HF_zx=HF_zx_sel(j,:);

        NQI_zz=NQI_zz_sel(j,:);
        NQI=zeros(Ni_ENDOR,3,3);

        NQI(:,:,:)=2*pi*NQI_sel(j,:,:,:);

        CS_zz=CS_zz_sel(j,:)*1e-12;


        S=S_sel(j);
        offsets=offsets_sel(j,:);

        endor_amp_tmp=zeros(1,Npts_EN);

        if abs(HF_zz(1))>1/spinSys("T2e")*0.1


        % loop over all offsets (spin manifolds)
        for i=1 : size(offsets,2)

            v_off_S=offsets(i);
            off_1=offsets(1);

            rho0=kehl_rho0(constants,paramsENDOR,B,geff,spinOps,spinSys,parameters,HF_zz,HF_zy,HF_zx,NQI_zz);

            %electron T1
            RT1e=kehl_relax_t1(rho0,Sx_D,spinSys("T1e"));
            RT1n=zeros(size(RT1e));
            RT2e=kehl_relax_t2(Sx_D,spinSys("T2e"));
            RT2n=zeros(size(RT2e));

            if N_spinSys==1
                for mm=1:Ni_ENDOR

                   %nuclear T1
                   RT1n=RT1e+kehl_relax_t1(rho0,Ix_D{mm},spinSys("T1n"));
                   RT2e=RT2e+kehl_relax_t2(Sx_D*Ix_D{mm},spinSys("T2dq"));
                   RT2n=RT2n+kehl_relax_t2(Ix_D{mm},spinSys("T2n"));
                end
            else
                RT1n=RT1e+kehl_relax_t1(rho0,Ix_D{1},spinSys("T1n"));
                RT2e=RT2e+kehl_relax_t2(Sx_D*Ix_D{1},spinSys("T2dq"));
                RT2n=RT2n+kehl_relax_t2(Ix_D{1},spinSys("T2n"));
            end

            R=RT1e+RT1n+RT2e+RT2n;

            % for tilted frame Rho relaxation times

            %electron T1
            RT1eRho=kehl_relax_t1(rho0,Sx_D,spinSys("T1eR"));
            RT1nRho=zeros(size(RT1eRho));

            % electron T2
            RT2eRho=kehl_relax_t2(Sx_D,spinSys("T2eR"));
            RT2nRho=zeros(size(RT2eRho));

            if N_spinSys==1
                for mm=1:Ni_ENDOR

                   %nuclear T1
                   RT1nRho=RT1eRho+kehl_relax_t1(rho0,Ix_D{mm},spinSys("T1nR"));

                   % double and zero quantum T2
                   RT2eRho=RT2eRho+kehl_relax_t2(Sx_D*Ix_D{mm},spinSys("T2dqR"));

                   % nuclear T2
                   RT2nRho=RT2nRho+kehl_relax_t2(Ix_D{mm},spinSys("T2nR"));
                end
            else

                %nuclear T1
                RT1nRho=RT1eRho+kehl_relax_t1(rho0,Ix_D{1},spinSys("T1nR"));

                % double and zero quantum T2
                RT2eRho=RT2eRho+kehl_relax_t2(Sx_D*Ix_D{1},spinSys("T2dqR"));

                % nuclear T2
                RT2nRho=RT2nRho+kehl_relax_t2(Ix_D{1},spinSys("T2nR"));
            end

            % full relaxation superoperator
            RRho=RT1eRho+RT1nRho+RT2eRho+RT2nRho;


            start_EN=paramsENDOR("start_EN");
            step_EN=paramsENDOR("step_EN");

            sl=expt("SL");
            oneE=expt("oneE");
            oneN=expt("oneN");

            Hfree_p=2*pi*v_off_S*Sz;
            if spinSys('N_SpinSys')>1
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
            else
                for mm=1:Ni_ENDOR
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
                end
            end

            % Integration step for the Signal to account for oscillation
            if v_off_S==0
                t9=1/(off_1*Nint);
            else
                t9=1/(v_off_S*Nint);
            end


            % loop over rf frequencies (x-axis)
            for a=1:Npts_EN
            % for a=1   %%% for testing

                % Radiofrequency
                v_RF=(start_EN+step_EN*(a-1));

                Hcorr=zeros(size(Hfree_p));
                HRF=Hfree_p;

                if parameters.Bterm==false
                    if spinSys('N_SpinSys')>1
                        m=1;
                        Hcorr=Hcorr+2*pi*v_RF*Iz{m};

                        HRF=HRF+2*pi*v_RF*Iz{m}+oneN*Iy{m};
                    else
                        for mm=1:Ni_ENDOR
                            Hcorr=Hcorr+2*pi*v_RF*Iz{mm};

                            HRF=HRF+2*pi*v_RF*Iz{mm}+oneN*Iy{mm};
                        end
                    end
                end

                Hfree=Hfree_p+Hcorr;

                HSL=Hfree+sl*Sy;

                [V,D]=eig(HSL);
                if ~issorted(diag(D))
                  [d,inds]=sort(diag(D));
                  D=diag(d);
                  V=V(:, inds);
                end
                %eigenvalues are order from negative large to positive large!
                W=V^(-1);

                HSL_t=W*HSL*W^(-1);


                HSL_t=full(hilb2liouv(sparse(HSL_t),'comm'));
                Hnonsel=full(hilb2liouv(sparse(Hfree+oneE*Sx),'comm'));

                Hfree=full(hilb2liouv(sparse(Hfree),'comm'));


                if parameters.Bterm==false

                    U3=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*full(hilb2liouv(sparse(HRF),'comm'))),t(3)));

                    U1=full(propagator(spinOps('spin_system'),1i*sparse(RRho-1i*HSL_t),t(1)));
                    U2=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hfree),t(2)));

                    U4=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hfree),t(4)));
                    U5=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hnonsel),t(5)));
                    U6=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hfree),t(6)));
                    U7=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hnonsel),t(7)));

                    U8=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hfree),t(6)+t(5)/2));
                    U9=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hfree),t9));
                else
                    U3=kehl_rf_bterm_rlx(parameters,expt,v_RF,Hfree_p,Iy,t(3),Ni_ENDOR,N_spinSys,R,spinOps('spin_system'));

                    U1=full(propagator(spinOps('spin_system'),1i*sparse(RRho-1i*HSL_t),t(1)));
                    U2=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hfree),t(2)));

                    U4=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hfree),t(4)));
                    U5=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hnonsel),t(5)));
                    U6=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hfree),t(6)));
                    U7=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hnonsel),t(7)));

                    U8=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hfree),t(6)+t(5)/2));
                    U9=full(propagator(spinOps('spin_system'),1i*sparse(R-1i*Hfree),t9));

                end

                % Evolve the densitymatrix

                rho=hilb2liouv(rho0,'statevec');

                rho_t=W*reshape(rho,sqrt(size(rho,1)),sqrt(size(rho,1)))*W^(-1);
                rho_t=hilb2liouv(rho_t,'statevec');
                rho_t=U1*rho_t;
                rho=V*reshape(rho_t,sqrt(size(rho_t,1)),sqrt(size(rho_t,1)))*V^(-1);
                rho=hilb2liouv(rho,'statevec');

                % rho = U1*rho;
                rho=U2*rho;
                rho=U3*rho;
                rho=U4*rho;
                rho=U5*rho;
                rho=U6*rho;
                rho=U7*rho;
                rho=U8*rho;


                value_Sy=0;
                for b=1:Nint
                   rho=U9*rho;
                   rho_f=reshape(rho,sqrt(size(rho,1)),sqrt(size(rho,1)));
                   value_Sy=value_Sy+real(trace(rho_f*Sy));
                end

                endor_amp_tmp(a)=endor_amp_tmp(a)+(value_Sy*S/(Nint*size(offsets,2)));

            end
        end
        end
        endor_amp=endor_amp+endor_amp_tmp;
    end
end

function kehl_spinlock_rlx_grumble(spin_system,parameters)
if (~isstruct(spin_system))||(~isfield(spin_system,'bas'))||(~isfield(spin_system,'comp'))
    error('spin_system must be a Spinach spin system structure.');
end
if ~isstruct(parameters)
    error('parameters must be a structure.');
end
end

function endor_amp=kehl_spinlock_calc(spin_system,parameters)

    % Check consistency
    kehl_spinlock_calc_grumble(spin_system,parameters);

    % Unpack legacy data
    constants=parameters.constants;
    spinOps=parameters.spinOps;
    spinSys=parameters.spinSys;
    expt=parameters.expt;
    paramsENDOR=parameters.paramsENDOR;
    EPR=parameters.epr;
    Sz=spinOps("Sz");
    Sx=spinOps("Sx");
    Sy=spinOps("Sy");
    Iz=spinOps("Iz");
    Ix=spinOps("Ix");
    Iy=spinOps("Iy");
    t=expt("t");
    Nint=8;

    N_spinSys=spinSys("N_SpinSys");

    Ni_ENDOR=spinSys("Ni_ENDOR");

    if N_spinSys>1
       nucs=Ni_ENDOR;
        for i=1:Ni_ENDOR-1
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

    I=spinSys("I");
    v_L=paramsENDOR("v_L");


    endor_amp=zeros(1,Npts_EN);
    if length(B_sel)==0
        % No resonance orientations were found
        return
    end


    % loop to repeat the calculation for every orientation
    parfor j=1:length(B_sel)

        % set parameters for this orientation
        geff=geff_sel(j);
        B=B_sel(j);

        HF_zz=HF_zz_sel(j,:);
        HF_zy=HF_zy_sel(j,:);
        HF_zx=HF_zx_sel(j,:);

        NQI_zz=NQI_zz_sel(j,:);

        NQI=zeros(Ni_ENDOR,3,3);

        NQI(:,:,:)=2*pi*NQI_sel(j,:,:,:);

        CS_zz=CS_zz_sel(j,:)*1e-12;
        D_zz=D_zz_sel(j,:);


        S=S_sel(j);
        offsets=offsets_sel(j,:);

        endor_amp_tmp=zeros(1,Npts_EN);


        % loop over all offsets (spin manifolds)
        for i=1 : size(offsets,2)

            v_off_S=offsets(i);
            off_1=offsets(1);

            rho0=kehl_rho0(constants,paramsENDOR,B,geff,spinOps,spinSys,parameters,HF_zz,HF_zy,HF_zx,NQI_zz);

            start_EN=paramsENDOR("start_EN");
            step_EN=paramsENDOR("step_EN");

            sl=expt("SL");
            oneE=expt("oneE");
            oneN=expt("oneN");

            Hfree_p=2*pi*v_off_S*Sz;
            if spinSys('N_SpinSys')>1
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

            else
                for mm=1:Ni_ENDOR
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
                end

                if spinSys("D_used")==true
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

            % mw pulses
            HSL_p=Hfree_p+sl*Sy;
            Hnonsel_p=Hfree_p+oneE*Sx;

            % Integration step for the Signal to account for oscillation
            if v_off_S==0
                t9=1/(off_1*Nint);
            else
                t9=1/(v_off_S*Nint);
            end


            % Calculate the propagators
            U1_p=full(propagator(spinOps('spin_system'),sparse(HSL_p),t(1)));
            U2_p=full(propagator(spinOps('spin_system'),sparse(Hfree_p),t(2)));

            U4_p=full(propagator(spinOps('spin_system'),sparse(Hfree_p),t(4)));
            U5_p=full(propagator(spinOps('spin_system'),sparse(Hnonsel_p),t(5)));
            U6_p=full(propagator(spinOps('spin_system'),sparse(Hfree_p),t(6)));
            U7_p=full(propagator(spinOps('spin_system'),sparse(Hnonsel_p),t(7)));
            U8_p=full(propagator(spinOps('spin_system'),sparse(Hfree_p),t(6)+t(5)/2));
            U9_p=full(propagator(spinOps('spin_system'),sparse(Hfree_p),t9));


            % loop over rf frequencies (x-axis)
            for a=1:Npts_EN
            % for a=1   %%% for testing

                % Radiofrequency
                v_RF=(start_EN+step_EN*(a-1));

                Hcorr=zeros(size(Hfree_p));
                HRF=Hfree_p;

                if parameters.Bterm==false
                    if spinSys('N_SpinSys')>1
                        m=1;
                        HRF=HRF+2*pi*v_RF*Iz{m}+oneN*Iy{m};
                    else
                        for mm=1:Ni_ENDOR
                            HRF=HRF+2*pi*v_RF*Iz{mm}+oneN*Iy{mm};
                        end
                    end
                end

                if spinSys('N_SpinSys')>1
                    m=1;
                    Hcorr=Hcorr+2*pi*v_RF*Iz{m};
                else

                    for mm=1:Ni_ENDOR
                        Hcorr=Hcorr+2*pi*v_RF*Iz{mm};
                    end
                end


                if parameters.Bterm==false
                    U3=full(propagator(spinOps('spin_system'),sparse(HRF),t(3)));

                    U1=U1_p*full(propagator(spinOps('spin_system'),sparse(Hcorr),t(1)));
                    U2=U2_p*full(propagator(spinOps('spin_system'),sparse(Hcorr),t(2)));

                    U4=U4_p*full(propagator(spinOps('spin_system'),sparse(Hcorr),t(4)));
                    U5=U5_p*full(propagator(spinOps('spin_system'),sparse(Hcorr),t(5)));
                    U6=U6_p*full(propagator(spinOps('spin_system'),sparse(Hcorr),t(6)));
                    U7=U7_p*full(propagator(spinOps('spin_system'),sparse(Hcorr),t(7)));

                    U8=U8_p*full(propagator(spinOps('spin_system'),sparse(Hcorr),t(6)+t(5)/2));
                    U9=U9_p*full(propagator(spinOps('spin_system'),sparse(Hcorr),t9));
                else
                    Hfree=Hfree_p;

                    U3=kehl_rf_bterm(parameters,expt,v_RF,Hfree,Iy,t(3),Ni_ENDOR,N_spinSys,spinOps('spin_system'));

                    U1=U1_p;
                    U2=U2_p;
                    U4=U4_p;
                    U5=U5_p;
                    U6=U6_p;
                    U7=U7_p;
                    U8=U8_p;
                    U9=U9_p;
                end

                HSL=HSL_p;
                if spinSys('N_SpinSys')>1
                    m=1;
                    HSL=HSL+2*pi*v_RF*Iz{m};
                else
                    for mm=1:Ni_ENDOR
                        HSL=HSL+2*pi*v_RF*Iz{mm};
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



                % Evolve the densitymatrix
                rho=rho0;
                rho=U1*rho*U1';
                rho_t=W*rho*W^(-1);
                rho_t=diag(diag(rho_t));
                rho=V*rho_t*V^(-1);
                rho=U2*rho*U2';
                rho=U3*rho*U3';

                if parameters.Bterm==true
                    rho=diag(diag(rho));
                end

                rho=U4*rho*U4';
                rho=U5*rho*U5';
                rho=U6*rho*U6';
                rho=U7*rho*U7';
                rho=U8*rho*U8';


                value_Sy=0;
                for b=1:Nint
                   rho=U9*rho*U9';
                   value_Sy=value_Sy+real(trace(rho*Sy));
                end

                if spinSys('N_SpinSys')>1
                    if parameters.powder==0
                        s=1;
                    else
                        s=size(offsets,2)/Ni_ENDOR;
                    end
                else
                    s=size(offsets,2);
                end

                endor_amp_tmp(a)=endor_amp_tmp(a)+(value_Sy*S/(Nint*s));

            end
        end
        endor_amp=endor_amp+endor_amp_tmp;
    end
end

function kehl_spinlock_calc_grumble(spin_system,parameters)
if (~isstruct(spin_system))||(~isfield(spin_system,'bas'))||(~isfield(spin_system,'comp'))
    error('spin_system must be a Spinach spin system structure.');
end
if ~isstruct(parameters)
    error('parameters must be a structure.');
end
end
