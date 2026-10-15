% Structure coefficient tables for the associative envelopes of su(mult)
% algebras. Syntax:
%
%    [product_table_left,product_table_right]=ist_product_table(mult)
%
% Parameters:
%
%    mult - multiplicity (number of energy levels) of
%           the spin in question
%
% Outputs:
%
%     product_table_left
%     product_table_right - structure coefficients in
%                           N-by-4 arrays of [n m k value] rows:
%
%           T{n}*T{m}=...+value*T{k}+...  (left table)
%
%           T{m}*T{n}=...+value*T{k}+...  (right table)
%
% Indices n, m, and k are one-based. These coefficients describe the IST
% expansion of left and right multiplicative action by T{n} on T{m}, as
% given in Eq 7.18 of the first edition of IK's book (normalisation is
% missing in the book, that's a typo). Numbering translation between
% single and double indices is given by lm2lin and lin2lm functions.
%
% Note: only nonzero coefficients are stored. Projection conservation,
%       the rank triangle, and odd rank sums at zero projections eliminate
%       analytically forbidden products before contraction. Coordinate
%       caches are separate from legacy dense caches; spin-1/2 is hard-coded.
%
% hannah.hogben@chem.ox.ac.uk
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=ist_product_table.m>

function [product_table_left,product_table_right]=ist_product_table(mult)

% Check consistency
grumble(mult);

% Hard-coded for spin-1/2 to avoid an expensive disk hit
if mult==2

    % Left side product table
    product_table_left=[1 1 1  1.0;
                        4 2 1 -0.5;
                        3 3 1  0.5;
                        2 4 1 -0.5;
                        2 1 2  0.5;
                        1 2 2  1.0;
                        3 2 2  0.5;
                        2 3 2 -0.5;
                        3 1 3  0.5;
                        4 2 3  0.5;
                        1 3 3  1.0;
                        2 4 3 -0.5;
                        4 1 4  0.5;
                        4 3 4  0.5;
                        1 4 4  1.0;
                        3 4 4 -0.5];

    % Right side product table
    product_table_right=[1 1 1  1.0;
                         4 2 1 -0.5;
                         3 3 1  0.5;
                         2 4 1 -0.5;
                         2 1 2  0.5;
                         1 2 2  1.0;
                         3 2 2 -0.5;
                         2 3 2  0.5;
                         3 1 3  0.5;
                         4 2 3 -0.5;
                         1 3 3  1.0;
                         2 4 3  0.5;
                         4 1 4  0.5;
                         4 3 4 -0.5;
                         1 4 4  1.0;
                         3 4 4  0.5];
    return;

end

% Generate cache record name
own_path=mfilename('fullpath');
own_path=own_path(1:(end-17));
table_file=[own_path 'ist_product_table_coords_' ...
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
    
    % Get tensor ranks and projections in the existing linear order
    [ranks,projs]=lin2lm(0:mult^2-1);

    % Bound storage by projection conservation before evaluating products
    capacity=0;
    for k=1:mult^2
        for m=1:mult^2
            capacity=capacity+nnz(projs==projs(k)-projs(m));
        end
    end
    product_table_left=zeros(capacity,4);
    product_table_right=zeros(capacity,4);
    left_count=0; right_count=0;

    % Populate only angular-momentum-allowed product coefficients
    for k=1:mult^2
        for m=1:mult^2
            for n=1:mult^2

                % Enforce projection conservation and the rank triangle
                if (projs(k)~=projs(n)+projs(m))||...
                   (ranks(k)<abs(ranks(n)-ranks(m)))||...
                   (ranks(k)>ranks(n)+ranks(m))
                    continue;
                end

                % Odd rank sums vanish when all three projections are zero
                if (projs(n)==0)&&(projs(m)==0)&&(projs(k)==0)&&...
                   (mod(ranks(n)+ranks(m)+ranks(k),2)~=0)
                    continue;
                end

                % Left product action: carefully tiptoe around extreme norms
                left_value=norms(n)*hdot((T{k}/norms(k)),...
                                        (T{n}/norms(n))*(T{m}/norms(m)));

                % Right product action: carefully tiptoe around extreme norms
                right_value=norms(n)*hdot((T{k}/norms(k)),...
                                         (T{m}/norms(m))*(T{n}/norms(n)));

                % Store nonzero coefficients without a numerical cutoff
                if left_value~=0
                    left_count=left_count+1;
                    product_table_left(left_count,:)=[n m k left_value];
                end
                if right_value~=0
                    right_count=right_count+1;
                    product_table_right(right_count,:)=[n m k right_value];
                end

            end
        end
    end

    % Discard unused capacity
    product_table_left=product_table_left(1:left_count,:);
    product_table_right=product_table_right(1:right_count,:);

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


