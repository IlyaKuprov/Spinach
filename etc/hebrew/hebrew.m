% IK's Hebrew flashcards function. The Excel file should contain
% one header row followed by three text columns:
%
%       A) English words
%       B) Hebrew translations
%       C) Latin transliterations
%
% It repeatedly shows a random Hebrew word, waits for <Enter>,
% then reveals the English word and its transliteration. Exit 
% with CTRL+C.
%
% ilya.kuprov@weizmann.ac.il

function hebrew(xlsFile)
    
    % Read the spreadsheet
    T=readtable(xlsFile,'TextType','string',      ...
                        'ReadVariableNames',true, ... 
                        'PreserveVariableNames',true);

    % Extract columns
    eng=T{:,1}; heb=T{:,2}; trans=T{:,3}; n=numel(eng);

    % Main flash card loop
    while true

        % Random row
        idx=randi(n);   

        % Show Hebrew
        disp(' '); disp(heb(idx));

        % Wait for user input, then reveal translation
        pause(); fprintf('%s (%s)\n',eng(idx),trans(idx));
        
    end

end

% "Oh don't worry, Ilya, we do have Wi-Fi
%  in the bomb shelter!"
%
% Michal Leskes to IK, in November 2023,
% as walls were shaking from missiles
% being intercepted over Rehovot

