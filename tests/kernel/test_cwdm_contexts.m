% Tests rejection of state-dependent chemistry in static contexts. Syntax:
%
%                    result=test_cwdm_contexts()
%
% Outputs:
%
%    result - callback-rate and bimolecular rejection checks for all
%             nine contexts that consume kinetics matrices
%
% Constant first-order kinetics must still reach a custom sequence. The
% sequence probes return the assembled matrix, without propagation.
%
% ilya.kuprov@weizmann.ac.il

function result=test_cwdm_contexts()

% Announce the static-context contract
fprintf('TESTING: CWDM static context chemistry boundaries\n');
result=new_test_result('kernel/cwdm_contexts','CWDM context boundaries',...
                      'Static contexts reject state-dependent reaction records.');
contexts={'liquid','imaging','crystal','powder','device',...
          'floquet','singlerot','doublerot','gridfree'};

% Exercise callback, mass-action, and constant kinetics at every boundary
for n=1:numel(contexts)
    for k=1:3

        % Construct two exchange pools and an inert third substance
        sys.magnet=1; sys.isotopes={'1H','1H','1H'};
        inter.chem.parts={1,2,3}; inter.chem.concs=[0.7 0.3 1];
        inter.chem.reactions={struct('reactants',1,'products',2,...
                                    'matching',[1 2],'rate',1)};
        if k==1, inter.chem.reactions{1}.rate=@(t)1+t; end
        if k==2
            inter.chem.reactions{1}.reactants=[1 3];
            inter.chem.reactions{1}.products=[2 3];
            inter.chem.reactions{1}.matching=[1 2;3 3];
        end
        bas.formalism='sphten-liouv'; bas.approximation={'none','none','none'};
        if strcmp(contexts{n},'device')
            sys.isotopes{3}='C3'; bas.formalism='zeeman-liouv';
        end
        s=test_spin_system(sys,inter,bas);

        % Supply valid small grids for each static context
        p.spins={'1H'}; p.offset=0; p.decouple={}; p.needs={};
        p.orientation=[0 0 0]; p.grid='single_crystal'; p.serial=true; p.verbose=0;
        p.rate=1000; p.axis=[0 0 1]; p.max_rank=1; p.tau_c=1e-3;
        p.rate_outer=800; p.rate_inner=2400; p.rank_outer=1; p.rank_inner=1;
        p.axis_outer=[0 0 1]; p.axis_inner=[0 0 1];
        p.npts=10; p.dims=0.01; p.deriv={'period',3}; p.diff=0; p.u=zeros(10,1);
        p.rlx_ph={}; p.rlx_op={}; p.rho0_ph={}; p.rho0_st={};
        p.coil_ph={}; p.coil_st={};
        if strcmp(contexts{n},'floquet'), p.grid='leb_2ang_rank_5'; end
        if strcmp(contexts{n},'gridfree'), p=rmfield(p,'grid'); end
        assumptions='nmr';
        if strcmp(contexts{n},'device'), assumptions='labframe'; end

        % Invoke production contexts and inspect the named boundary
        rejected=false; reached=false;
        try
            if strcmp(contexts{n},'imaging')
                answer=imaging(s,@spatial_probe,p);
            else
                answer=feval(contexts{n},s,@matrix_probe,p,assumptions);
            end
            reached=isfinite(norm(answer,1));
        catch err
            rejected=strcmp(err.identifier,['Spinach:' contexts{n} ':stateDependentKinetics'])&&...
                     contains(err.message,'state-dependent reaction records')&&...
                     contains(err.message,'step/iserstep')&&...
                     contains(err.message,'examples/kinetics/nonlinear/bimolecular_closures.m')&&...
                     contains(err.message,'examples/microfluidics/reacting_flow_nmr.m');
            if ~rejected, fprintf('CONTEXT_UNEXPECTED %s %s\n',contexts{n},err.message); end
        end
        if k<3
            result=test_true(result,[contexts{n} ' rejection ' num2str(k)],rejected,...
                             'state-dependent kinetics must receive the actionable context error');
        else
            result=test_true(result,[contexts{n} ' constant kinetics'],reached,...
                             'constant matrix kinetics must still reach the pulse sequence');
        end
    end
end

end

% Return the static spin or rotor generator for the acceptance control
function answer=matrix_probe(spin_system,parameters,H,R,K) %#ok<INUSD>
answer=full(K);
end

% Return the static imaging generator for the acceptance control
function answer=spatial_probe(spin_system,parameters,H,R,K,G,F) %#ok<INUSD>
answer=full(K);
end


