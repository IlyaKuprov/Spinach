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
% rallelisation call stack.
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
                                        F(job_number)=parfeval(@brw_compute_kernel,0,spin_system,...
                                                               weights{s}(j),job_number,upper_limit); %#ok<SAGROW>
                                                               
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
    
% Intermediate clean-up
clear('Q','L0','A','B','C','D','weights','rates','states');

% Gather data from the nodes
if exist('F','var')
    
    % Wait for the asynchronous calculation to finish
    report(spin_system,'waiting for data from the nodes...'); 
    drawnow(); wait(F);
    
    % Retrieve from ValueStore
    store_keys=cell(numel(F),1);
    for n=1:numel(F)

        % Make ValueStore keys
        if ~isempty(F(n).Error)

            % Rethrow remote error
            rethrow(F(n).Error.remotecause{1});

        else

            % Build the ValueStore key
            store_keys{n}=['redfield_int_batch_' num2str(n)];

        end

    end

    % Retrieve from ValueStore, assemble, and add to the total
    XYZ=cell2mat(get(store,store_keys)); 
    remove(store,store_keys); dim=size(R,1);
    R=R+sparse(XYZ(:,1),XYZ(:,2),XYZ(:,3),dim,dim); 
    report(spin_system,'Redfield superoperator done.');

    % Final clean-up
    clear('F','XYZ');

else

    % A rare situation when there is nothing significant in the result
    report(spin_system,'WARNING - Redfield module found no significant elements in Q.');

end

% BRW integral kernel for asynchronous execution
function brw_compute_kernel(spin_system,w,job_number,upper_lim)

% Get the data from the ValueStore
store=getCurrentValueStore();
store_keys={['brw_integrator_batch_' num2str(job_number) '_A'],...
            ['brw_integrator_batch_' num2str(job_number) '_B'],...
            ['brw_integrator_batch_' num2str(job_number) '_C'],...
            ['brw_integrator_batch_' num2str(job_number) '_D']};
ABCD=get(store,store_keys); remove(store,store_keys);

% Compute the integral and do an intermediate clean-up
R=-w*ABCD{1}*expmint(spin_system,ABCD{2},ABCD{3},ABCD{4},upper_lim); 
clear('ABCD'); R=clean_up(spin_system,R,spin_system.tols.rlx_zero);

% Convert CSR format into XYZ for adding
[row,col,val]=find(R); R=[row col val];

% Send to the pool ValueStore
put(store,{['redfield_int_batch_' num2str(job_number)]},{R}); 

end

% Owning a business means working 80 hours a week 
% so you can avoid working 40 hours for someone else.
%
% Ramona Arnett

