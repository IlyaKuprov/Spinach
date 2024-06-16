% Redfield integral evaluation, the asynchronous parallel path. This
% include is called from within relaxation.m and follows the notati-
% on used in IK's paper:
%
%           http://dx.doi.org/10.1016/j.jmr.2010.12.004
%
% with the difference that the numerical quadrature method proposed
% there has been superceded by the much faster auxiliary matrix me-
% thod described in:
%
%                http://dx.doi.org/10.1063/1.4928978
%
% This include is called when relaxation.m is at the top of the pa-
% rallelisation call stack. When relaxation.m is called from inside 
% a parallel loop, the serial version of this include is used.
%
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=redfield_integral_async.m>

% Async job counter
job_number=0;

% Sum over spherical ranks
for n=1:numel(Q)
    
    % Sum over first projection index
    for k=1:(2*n+1)
        
        % Sum over second projection index
        for m=1:(2*n+1)
            
            % Only proceed if necessary
            if cheap_norm(Q{n}{k,m})>0
                
                % Sum over third projection index
                for p=1:(2*n+1)
                    
                    % Sum over fourth projection index
                    for q=1:(2*n+1)
                        
                        % Only proceed if necessary
                        if cheap_norm(Q{n}{p,q})>0
                            
                            % Get decay weights and rates for correlation function
                            [weights,rates,states]=corrfun(spin_system,n,k,m,p,q);
                            
                            % Loop over chemical species
                            for s=1:numel(states)
                                
                                % Loop over correlation function exponentials
                                for j=1:numel(weights{s})
                                    
                                    % Only proceed if necessary
                                    if abs(weights{s}(j))>0

                                        % Increment the job counter
                                        job_number=job_number+1;
                                        
                                        % Set the upper integration limit according to the accuracy goal
                                        upper_limit=-1.5*(1/rates{s}(j))*log(1/spin_system.tols.rlx_integration);
                    
                                        % Inform the user
                                        report(spin_system,['queuing up L=' num2str(n)  ', k=' num2str(n+1-k,'%+d') ... 
                                                            ', m=' num2str(n+1-m,'%+d') ', p=' num2str(n+1-p,'%+d') ...
                                                            ', q=' num2str(n+1-q,'%+d') ', chemical species ' num2str(s) '...']);
                                            
                                        % Kill the terms in L0 that are irrelevant on the time scale of the integration
                                        B=clean_up(spin_system,L0,spin_system.tols.rlx_integration/abs(upper_limit));
                                                            
                                        % Prepare the relevant matrices
                                        A=Q{n}{k,m}; C=Q{n}{p,q}'; D=B-1i*rates{s}(j)*speye(size(B));
                                        
                                        % Obliterate irrelevant elements
                                        A(~states{s},~states{s})=0; B(~states{s},~states{s})=0;
                                        C(~states{s},~states{s})=0; D(~states{s},~states{s})=0;
                                            
                                        % Put the matrices on the ValueStore
                                        store=gcp('nocreate').ValueStore;
                                        store.put({['brw_integrator_batch_' num2str(job_number) '_A'],...
                                                   ['brw_integrator_batch_' num2str(job_number) '_B'],...
                                                   ['brw_integrator_batch_' num2str(job_number) '_C'],...
                                                   ['brw_integrator_batch_' num2str(job_number) '_D']},{A,B,C,D});
                                       
                                        % Queue up Redfield integral using the auxiliary matrix method
                                        F(job_number)=parfeval(@brw_compute_kernel,1,spin_system,...
                                                               weights{s}(j),job_number,upper_limit); %#ok<SAGROW>
                                                               
                                    end
                                end
                            end
                        end             
                    end                 % Everything that needs to be said has        
                end                     % already been said. But, since no one        
            end                         % was listening, everything must be        
        end                             % said again. - Andre Gide
    end
end
    
% Intermediate clean-up
clear('Q','L0','A','B','C','D','weights','rates','states');

% Gather data from the nodes
if exist('F','var')

    % Tell the user we are retrieving the data now
    report(spin_system,'getting data from workers...'); 
    drawnow(); XYZ=cell(job_number,1);

    % Retrieve the data
    for n=1:job_number

        % Get the result ID
        job_id=fetchNext(F);

        % Make ValueStore key
        if ~isempty(F(job_id).Error)

            % Rethrow remote error if we bombed out
            rethrow(F(job_id).Error.remotecause{1});

        else

            % Build ValueStore key for the job output
            store_key=['redfield_int_batch_' num2str(job_id)];

        end

        % Retrieve the fragment and delete it from the ValueStore
        XYZ{n}=cell2mat(get(store,store_key)); remove(store,store_key); 
        
    end

    % Assemble and add to the total
    clear('F'); XYZ=cell2mat(XYZ); dim=size(R,1);
    R=R+sparse(XYZ(:,1),XYZ(:,2),XYZ(:,3),dim,dim); 

    % Final clean-up
    clear('XYZ','store');

else

    % A rare situation when there is nothing significant in the result
    report(spin_system,'WARNING - Redfield module found no significant elements in Q.');

end

% Asynchronous BRW integral kernel, careful memory management
function job_id=brw_compute_kernel(spin_system,w,job_id,upper_lim)

% Get data from pool ValueStore
store=getCurrentValueStore();
store_keys={['brw_integrator_batch_' num2str(job_id) '_A'],...
            ['brw_integrator_batch_' num2str(job_id) '_B'],...
            ['brw_integrator_batch_' num2str(job_id) '_C'],...
            ['brw_integrator_batch_' num2str(job_id) '_D']};
ABCD=get(store,store_keys); remove(store,store_keys);

% Compute the integral and do an intermediate clean-up
R=-w*ABCD{1}*expmint(spin_system,ABCD{2},ABCD{3},ABCD{4},upper_lim); 
clear('ABCD'); R=clean_up(spin_system,R,1e-2*spin_system.tols.rlx_zero);

% Convert CSR sparse format into XYZ format
[row,col,val]=find(R); clear('spin_system','R');

% Send to the pool ValueStore and delete local copy
put(store,{['redfield_int_batch_' num2str(job_id)]},...
          {[row col val]}); 
clear('row','col','val','store');

end

% Owning a business means working 80 hours a week 
% so you can avoid working 40 hours for someone else.
%
% Ramona Arnett

