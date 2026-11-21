% Test: kernel/cce_erbium_s12 | Er:CaWO4 secular pair reduction
% Compare the analytical Er:CaWO4 CCE-2 pair product to dense, full pair
% Hahn propagation. Returns the standard Spinach test result record.
%
% ilya.kuprov@weizmann.ac.il

function result=test_cce_erbium_s12()

% Three tungsten sites in the region of resolvable flip-flop decay
result=new_test_result('kernel/cce_erbium_s12','Er:CaWO4 secular pair reduction',...
                       'Compare the analytical pair product with explicit Hahn propagation.');
xyz=1e-9*[4.456465 0.000000 -0.2843425;...
          4.194320 0.262145  0.2843425;...
          4.456465 0.262145  0.0000000];
field_axis=[cosd(46.5);sind(46.5);0];
gamma_w=spin('183W');
hyperfine=1e-7*8.38*9.2740100657e-24*gamma_w;
dipolar=1e-7*(6.62607015e-34/(2*pi))*gamma_w^2;
times=linspace(0,0.060,61).';

% Construct full four-dimensional secular pair Hamiltonians
sx=[0 1;1 0]/2;
sy=[0 -1i;1i 0]/2;
sz=diag([1 -1])/2;
pair_op=kron(sz,sz)-(kron(sx,sx)+kron(sy,sy))/2;
zeeman=-gamma_w*0.067*(kron(sz,eye(2))+kron(eye(2),sz));
reference=ones(size(times));
for j=1:size(xyz,1)-1
    for k=j+1:size(xyz,1)

        % Compute dipolar tensor projections independently of the analytical function
        r=xyz(k,:)-xyz(j,:);
        u=r.'/norm(r);
        D=dipolar*(eye(3)-3*(u*u.'))/norm(r)^3;
        r=xyz(j,:);
        u=r.'/norm(r);
        A=hyperfine*(eye(3)-3*(u*u.'))/norm(r)^3;
        r=xyz(k,:);
        u=r.'/norm(r);
        B=hyperfine*(eye(3)-3*(u*u.'))/norm(r)^3;
        coupling=(field_axis.'*A*field_axis)*kron(sz,eye(2))+...
                 (field_axis.'*B*field_axis)*kron(eye(2),sz);
        H=zeeman+(field_axis.'*D*field_axis)*pair_op;
        for n=1:numel(times)

            % Compute the trace of the full Hahn commutator including spectator levels
            ua=expm(-1i*(H+coupling/2)*times(n)/2);
            ub=expm(-1i*(H-coupling/2)*times(n)/2);
            reference(n)=reference(n)*trace((ua*ub)'*(ub*ua))/4;

        end

    end
end

% Compare the analytical product to the explicit propagators
answer=cce_secular_pairs(xyz,field_axis,hyperfine,dipolar,times,0);
result=test_true(result,'Resolved decay',min(real(reference))<0.999,...
                 'The finite bath must distinguish the echo from unity.');
result=test_close(result,'Secular pair product',answer,reference,1e-8,1e-8,...
                  'Analytical flip-flop blocks must match dense propagation.');

% A singleton and a commuting nuclear bath both have perfect Hahn refocusing
for k=[1 size(xyz,1)]

    % Vanishing detuning exercises rows with no retained factors
    answer=cce_secular_pairs(xyz(1:k,:),field_axis,0,dipolar,times,0);
    result=test_true(result,sprintf('Commuting case (%d spins)',k),all(answer==1),...
                     'Without hyperfine detuning the Hahn echo must equal unity.');

end

end


