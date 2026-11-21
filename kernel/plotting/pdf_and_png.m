% Drops publisher quality PDF and PNG files of the current
% figure into the current directory. Syntax:
%
%                  pdf_and_png(file_name)
%
% Parameters:
%
%    file_name - character string without extension
%
% Ouput:
%
%    this function writes two files
%
% ilya.kuprov@weizmann.ac.il

function pdf_and_png(file_name)

% Write the PDF file
exportgraphics(gcf,[file_name '.pdf'],'ContentType','Vector');

% Write the PNG file
exportgraphics(gcf,[file_name '.png'],'Resolution',1200);

end

% If it doesn't work horizontally as prose,
% it probably
% won't work
% any better
% vertically
% pretending
% to be poetry.
%
% Robert Brault

