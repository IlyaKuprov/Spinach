% IK's Hebrew flashcards function. The Excel files should contain
% Hebrew vocabulary in two separate spreadsheets:
%
%       nouns.xlsx - English, masculine singular, feminine singular,
%                    masculine plural, feminine plural, invariant,
%                    notes
%       verbs.xlsx - English, infinitive, masculine singular,
%                    feminine singular, masculine plural,
%                    feminine plural, notes
%
% Syntax:
%
%       hebrew()
%       hebrew(mode)
%       ncards=hebrew(mode,max_cards)
%
% The mode may be 'forward', 'backward', or 'both'. In forward mode,
% English is shown first and Hebrew is revealed after <Enter>. In
% backward mode, Hebrew is shown first and English is revealed after
% <Enter>. In both mode, the direction is random on every card. Exit
% an open-ended run with CTRL+C.
%
% ilya.kuprov@weizmann.ac.il

function ncards=hebrew(mode,max_cards) %#NHEAD

    % Set default direction mode
    if nargin<1
        mode="both";
    end

    % Set default number of cards
    if nargin<2
        max_cards=inf;
    end

    % Normalise and validate inputs
    mode=lower(string(mode));
    grumble(mode,max_cards);

    % Locate the spreadsheet directory
    this_file=mfilename('fullpath');
    root_dir=fileparts(this_file);

    % Load all usable flashcards
    cards=load_cards(root_dir);
    ncards=height(cards);

    % Return immediately for spreadsheet-loading tests
    if max_cards==0
        return
    end

    % Start with a random card order
    card_order=randperm(ncards);
    order_pos=1; shown=0;

    % Main flash card loop
    while shown<max_cards

        % Reshuffle after all cards have been used
        if order_pos>ncards
            card_order=randperm(ncards);
            order_pos=1;
            fprintf('\nAll cards shown; reshuffling.\n');
        end

        % Select the next card
        card_idx=card_order(order_pos);
        order_pos=order_pos+1;
        shown=shown+1;

        % Randomise direction in mixed mode
        card_mode=mode;
        if mode=="both"
            if rand()<0.5
                card_mode="forward";
            else
                card_mode="backward";
            end
        end

        % Show the selected card
        show_card(cards(card_idx,:),card_mode);

    end

end


function cards=load_cards(root_dir)

    % Start from an empty table
    cards=table(strings(0,1),strings(0,1),strings(0,1),strings(0,1),...
                'VariableNames',{'english','hebrew','form','source'});

    % Load noun/adjective forms
    noun_file=fullfile(root_dir,'nouns.xlsx');
    noun_forms=["masculine singular","feminine singular",...
                "masculine plural","feminine plural","invariant"];
    cards=[cards;read_cards(noun_file,"noun",noun_forms,2:6)];

    % Load verb forms
    verb_file=fullfile(root_dir,'verbs.xlsx');
    verb_forms=["infinitive","masculine singular","feminine singular",...
                "masculine plural","feminine plural"];
    cards=[cards;read_cards(verb_file,"verb",verb_forms,2:6)];

    % Refuse to run on empty input
    if isempty(cards)
        error('No flashcards found in nouns.xlsx and verbs.xlsx.');
    end

end


function cards=read_cards(file_name,source,form_names,form_cols)

    % Read the spreadsheet as strings
    T=readtable(file_name,'TextType','string','PreserveVariableNames',true);

    % Extract and clean English prompts
    english=clean_text(T{:,1});

    % Initialise card arrays
    eng_cards=strings(0,1); heb_cards=strings(0,1);
    form_cards=strings(0,1); src_cards=strings(0,1);

    % Loop over grammatical form columns
    for n=1:numel(form_cols)

        % Extract the current Hebrew form
        hebrew=clean_text(T{:,form_cols(n)});

        % Loop over spreadsheet rows
        for k=1:numel(english)

            % Keep only complete Hebrew flashcard entries
            if strlength(english(k))>0&&strlength(hebrew(k))>0&&contains_hebrew(hebrew(k))
                eng_cards(end+1,1)=english(k); %#ok<AGROW>
                heb_cards(end+1,1)=hebrew(k); %#ok<AGROW>
                form_cards(end+1,1)=form_names(n); %#ok<AGROW>
                src_cards(end+1,1)=source; %#ok<AGROW>
            end

        end

    end

    % Build the card table
    cards=table(eng_cards,heb_cards,form_cards,src_cards,...
                'VariableNames',{'english','hebrew','form','source'});

end


function text=clean_text(text)

    % Convert spreadsheet cells into trimmed strings
    text=string(text);
    text(ismissing(text))="";
    text=strtrim(text);

end


function tf=contains_hebrew(text)

    % Detect Hebrew code points
    tf=~isempty(regexp(char(text),'[\x{0590}-\x{05FF}]','once'));

end


function show_card(card,mode)

    % Build the English side of the card
    eng_text=sprintf('%s [%s, %s]',char(card.english),...
                     char(card.source),char(card.form));

    % Show the prompt and wait for Enter
    fprintf('\n');
    if mode=="forward"
        fprintf('%s\n',eng_text);
        input('', 's');
        fprintf('%s\n',char(card.hebrew));
    else
        fprintf('%s\n',char(card.hebrew));
        input('', 's');
        fprintf('%s\n',eng_text);
    end

end


function grumble(mode,max_cards)
if ~isscalar(mode)||~ismember(mode,["forward","backward","both"])
    error('mode must be ''forward'', ''backward'', or ''both''.');
end
if ~isnumeric(max_cards)||~isscalar(max_cards)||max_cards<0
    error('max_cards must be a non-negative real scalar.');
end
end

% "Oh don't worry, Ilya, we do have Wi-Fi
%  in the bomb shelter!"
%
% Michal Leskes to IK, in November 2023,
% as walls were shaking from missiles
% being intercepted over Rehovot

% #NWIKI


