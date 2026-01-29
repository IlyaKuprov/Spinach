% Structure coefficient tables for single transition operators. Syntax:
%
%              [pt_left,pt_right]=st_product_table(nlevels)
%
% The input parameter is the dimension of the density matrix. The out-
% put contains structure coefficients in the following conventions
% (S{m} and S{k} are normalised):
%
%                 S{n}*S{m}=...+pt_left(n,m,k)*S{k}+...
%
%                 S{m}*S{n}=...+pt_right(n,m,k)*S{k}+...
%
% corresponding to the expansion of the left and the right multiplica-
% tive action by S{n} on S{m} as given in Eq 7.18 of the first edition
% of IK's book (normalisation is missing in the book, that's a typo).
% Numbering translation between single and double index is given by 
% kq2lin and lin2kq functions.
%
% Note: these are expensive tables, disk cache is used automatically.
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=st_product_table.m>

function [pt_left,pt_right]=st_product_table(nlevels)

% Check consistency
grumble(nlevels);

% Generate cache record name
own_path=mfilename('fullpath');
own_path=own_path(1:(end-16));
table_file=[own_path 'st_product_table_' ...
            num2str(nlevels) '.mat'];

% Check the cache
if exist(table_file,'file')
    
    % Lift data from the cache if the file is already available
    load(table_file,'pt_left','pt_right');
    
else
    
    % Get ST operators
    B=sin_tran(nlevels);

    % Preallocate the arrays
    pt_left= zeros(nlevels^2,nlevels^2,nlevels^2);
    pt_right=zeros(nlevels^2,nlevels^2,nlevels^2);
    
    % Get the structure coefficients
    for k=1:nlevels^2
        for m=1:nlevels^2
            for n=1:nlevels^2

                % Left product action
                pt_left(n,m,k)= hdot(B{k},B{n}*B{m});

                % Right product action
                pt_right(n,m,k)=hdot(B{k},B{m}*B{n});

            end
        end
    end

    try % Try to save a cache record, but don't insist
        save(table_file,'pt_left','pt_right'); drawnow;
    catch
        warning('Spinach installation appears to be write-protected');
    end

end

end

% Consistency enforcement
function grumble(nlevels)
if (~isnumeric(nlevels))||(~isreal(nlevels))||...
   (~isscalar(nlevels))||(nlevels<1)||(mod(nlevels,1)~=0)
    error('mult must be a positive real integer.');
end
end

% Из толстых книг нельзя узнать ничего нового. Толстые 
% книги - это кладбище, в котором погребены отслужившие
% свой век идеи прошлого.
%
% Лев Ландау

