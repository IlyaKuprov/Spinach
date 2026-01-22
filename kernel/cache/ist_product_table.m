% Structure coefficient tables for the associative envelopes of su(mult)
% algebras. Syntax:
%
%    [product_table_left,product_table_right]=ist_product_table(mult)
%
% The input parameter is the multiplicity of the spin in question. Output
% contains the structure coefficients in the following conventions (T{m}
% and T{k} are normalised):
%
%            T{n}*T{m}=...+product_table_left(n,m,k)*T{k}+...
%
%            T{m}*T{n}=...+product_table_right(n,m,k)*T{k}+...
%
% corresponding to the IST expansion of the left and the right multiplica-
% tive action by T{n} on T{m} as given in Eq 7.18 of the first edition of
% IK's book (normalisation is missing in the book, that's a typo). Number-
% ing translation between single and double index is given by lm2lin and
% lin2lm functions.
%
% Note: these are expensive tables, disk cache is used automatically,
%       the smallest and most frequently used case is hard-coded.
%
% hannah.hogben@chem.ox.ac.uk
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=ist_product_table.m>

function [product_table_left,product_table_right]=ist_product_table(mult)

% Check consistency
grumble(mult);

% Hard-coded for spin-1/2 to 
% avoid expensive disk hit
if mult==2

    % Left side product table
    product_table_left=zeros(4,4,4);
    product_table_left(:,:,1)= [1.0     0     0     0
                                  0     0     0  -0.5
                                  0     0   0.5     0
                                  0  -0.5     0     0];
    product_table_left(:,:,2)= [  0   1.0     0     0
                                0.5     0  -0.5     0
                                  0   0.5     0     0
                                  0     0     0     0];
    product_table_left(:,:,3)= [  0     0   1.0     0
                                  0     0     0  -0.5
                                0.5     0     0     0
                                  0   0.5     0     0];
    product_table_left(:,:,4)= [  0     0     0   1.0
                                  0     0     0     0
                                  0     0     0  -0.5
                                0.5     0   0.5     0];

    % Right side product table
    product_table_right=zeros(4,4,4);
    product_table_right(:,:,1)=[1.0     0     0     0
                                  0     0     0  -0.5
                                  0     0   0.5     0
                                  0  -0.5     0     0];
    product_table_right(:,:,2)=[  0   1.0     0     0
                                0.5     0   0.5     0
                                  0  -0.5     0     0
                                  0     0     0     0];
    product_table_right(:,:,3)=[  0     0   1.0     0
                                  0     0     0   0.5
                                0.5     0     0     0
                                  0  -0.5     0     0];
    product_table_right(:,:,4)=[  0     0     0   1.0
                                  0     0     0     0
                                  0     0     0   0.5
                                0.5     0  -0.5     0]; return; 

end

% Generate cache record name
own_path=mfilename('fullpath');
own_path=own_path(1:(end-17));
table_file=[own_path 'ist_product_table_' ...
            num2str(mult) '.mat'];

% Check the cache
if exist(table_file,'file')
    
    % Lift data from the cache if the file is already available
    load(table_file,'product_table_left','product_table_right');
    
else
    
    % Get IST operators
    T=irr_sph_ten(mult);

    % Precompute norms
    norms=zeros(mult^2,1);
    for n=1:mult^2
        norms(n)=sqrt(hdot(T{n},T{n}));
    end
    
    % Preallocate IST product tables
    product_table_left= zeros(mult^2,mult^2,mult^2);
    product_table_right=zeros(mult^2,mult^2,mult^2);

    % Populate product tables
    for k=1:mult^2
        for m=1:mult^2
            for n=1:mult^2

                % Left product action: carefully tiptoe around extreme norms
                product_table_left(n,m,k)= norms(n)*hdot((T{k}/norms(k)),...
                                                        (T{n}/norms(n))*...
                                                        (T{m}/norms(m)));

                % Right product action: carefully tiptoe around extreme norms
                product_table_right(n,m,k)=norms(n)*hdot((T{k}/norms(k)),...
                                                        (T{m}/norms(m))*...
                                                        (T{n}/norms(n)));

            end
        end
    end

    try % Try to save a cache record, but do not insist
        save(table_file,'product_table_left',...
                        'product_table_right','-v7.3'); drawnow;
    catch
        warning('Spinach directory appears to be write-protected');
    end

end

end

% Consistency enforcement
function grumble(mult)
if (numel(mult)~=1)||(~isnumeric(mult))||...
   (~isreal(mult))||(mult<2)||(mod(mult,1)~=0)
    error('mult must be a real integer greater than 1.');
end
end

% According to Oxford Chemistry folklore, Peter Atkins has once asked the
% following question at an interview for a Lecturer post: "What is it that
% you have done that a technician would not do?". The candidate produced a
% reasonable answer to the effect that a technician would not have the re-
% quired skills. A better answer was suggested by an interview committee 
% member some time later: "And what is it that you, Atkins, have done that
% a journalist would not do?"

