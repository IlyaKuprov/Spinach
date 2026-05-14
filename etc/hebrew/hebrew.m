% IK's Hebrew flashcards function. The Excel files should contain
% Hebrew vocabulary in separate spreadsheets:
%
%       nouns.xlsx - English, masculine singular, feminine singular,
%                    masculine plural, feminine plural, invariant,
%                    notes
%       adverbs.xlsx, directions.xlsx, greetings.xlsx, languages.xlsx,
%                    numbers.xlsx, particles.xlsx, phrases.xlsx,
%                    prepositions.xlsx, pronouns.xlsx,
%                    proper_nouns.xlsx, quantifiers.xlsx,
%                    sex_terms.xlsx, weekdays.xlsx - same columns as
%                    nouns.xlsx
%       adjectives.xlsx - English, masculine singular, feminine singular,
%                         masculine plural, feminine plural, notes
%       question_words.xlsx - English, masculine singular,
%                             feminine singular, masculine plural,
%                             feminine plural, invariant, notes
%       verbs.xlsx - English, infinitive, masculine singular,
%                    feminine singular, masculine plural,
%                    feminine plural, notes
%
% Syntax:
%
%       hebrew()
%       hebrew(mode)
%       ncards=hebrew(mode,max_cards)
%       hebrew('gui')
%
% The mode may be 'forward', 'backward', 'both', or 'gui'. In forward mode,
% English is shown first and Hebrew is revealed after <Enter>. In
% backward mode, Hebrew is shown first and English is revealed after
% <Enter>. In both mode, the direction is random on every card. Exit
% an open-ended run with CTRL+C. In gui mode, a one-button graphical
% flashcard window is opened; the button reveals the answer on the
% first click and advances to another randomly selected card on the
% second click.
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

    % Start the graphical flashcard interface
    if mode=="gui"
        show_gui(cards);
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

    % Define the forms used by noun-like word classes
    noun_forms=["masculine singular","feminine singular",...
                "masculine plural","feminine plural","invariant"];

    % List noun-like word class files
    lex_files={...
        'nouns.xlsx',"noun";...
        'proper_nouns.xlsx',"proper noun";...
        'languages.xlsx',"language";...
        'weekdays.xlsx',"weekday";...
        'numbers.xlsx',"number";...
        'pronouns.xlsx',"pronoun";...
        'prepositions.xlsx',"preposition";...
        'adverbs.xlsx',"adverb";...
        'particles.xlsx',"particle";...
        'quantifiers.xlsx',"quantifier";...
        'directions.xlsx',"direction";...
        'phrases.xlsx',"phrase";...
        'greetings.xlsx',"greeting";...
        'sex_terms.xlsx',"sex term";...
        'question_words.xlsx',"question word"};

    % Load noun-like and invariant word classes
    for n=1:size(lex_files,1)
        lex_file=fullfile(root_dir,lex_files{n,1});
        cards=[cards;read_cards(lex_file,lex_files{n,2},noun_forms,2:6)];
    end

    % Load adjective forms
    adj_file=fullfile(root_dir,'adjectives.xlsx');
    adj_forms=["masculine singular","feminine singular",...
               "masculine plural","feminine plural"];
    cards=[cards;read_cards(adj_file,"adjective",adj_forms,2:5)];

    % Load verb forms
    verb_file=fullfile(root_dir,'verbs.xlsx');
    verb_forms=["infinitive","masculine singular","feminine singular",...
                "masculine plural","feminine plural"];
    cards=[cards;read_cards(verb_file,"verb",verb_forms,2:6)];

    % Refuse to run on empty input
    if isempty(cards)
        error('No flashcards found in the Hebrew vocabulary spreadsheets.');
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


function show_gui(cards)

    % Select a Hebrew-capable display font
    gui_font=select_font();

    % Create the flashcard window
    gui_fig=figure('Name','Hebrew flashcards','NumberTitle','off',...
                   'MenuBar','none','ToolBar','none','Color',[1 1 1],...
                   'Units','pixels','Position',[100 100 760 430]);

    % Create the card direction label
    context_ctrl=uicontrol(gui_fig,'Style','text','Units','normalized',...
                           'Position',[0.05 0.88 0.90 0.06],...
                           'BackgroundColor',[1 1 1],'FontName',char(gui_font),...
                           'FontSize',12,'HorizontalAlignment','center');

    % Create the current-word label
    prompt_ctrl=uicontrol(gui_fig,'Style','text','Units','normalized',...
                          'Position',[0.05 0.54 0.90 0.30],...
                          'BackgroundColor',[1 1 1],'FontName',char(gui_font),...
                          'FontSize',34,'FontWeight','bold',...
                          'HorizontalAlignment','center');

    % Create the answer label
    answer_ctrl=uicontrol(gui_fig,'Style','text','Units','normalized',...
                          'Position',[0.05 0.29 0.90 0.20],...
                          'BackgroundColor',[1 1 1],'FontName',char(gui_font),...
                          'FontSize',24,'HorizontalAlignment','center');

    % Create the single control button
    button_ctrl=uicontrol(gui_fig,'Style','pushbutton','Units','normalized',...
                          'Position',[0.35 0.08 0.30 0.13],...
                          'FontName',char(gui_font),'FontSize',16,...
                          'Callback',@button_press);

    % Initialise the GUI state
    gui_state.cards=cards;
    gui_state.card_idx=0;
    gui_state.revealed=false;
    gui_state.prompt="";
    gui_state.answer="";
    gui_state.context="";

    % Display the first random card
    next_card();

    function button_press(~,~)

        % Reveal or advance depending on the current state
        if gui_state.revealed
            next_card();
        else
            gui_state.revealed=true;
            set(answer_ctrl,'String',char(gui_state.answer));
            set(button_ctrl,'String','Next');
        end

    end

    function next_card()

        % Select a random card without immediate repetition where possible
        new_idx=randi(height(gui_state.cards));
        if height(gui_state.cards)>1&&new_idx==gui_state.card_idx
            new_idx=mod(new_idx,height(gui_state.cards))+1;
        end

        % Randomise the question direction
        if rand()<0.5
            card_mode="forward";
        else
            card_mode="backward";
        end

        % Prepare the prompt and answer strings
        gui_state.card_idx=new_idx;
        gui_state.revealed=false;
        [gui_state.prompt,gui_state.answer,gui_state.context]=...
            card_sides(gui_state.cards(new_idx,:),card_mode);

        % Refresh the visible card controls
        set(context_ctrl,'String',char(gui_state.context));
        set(prompt_ctrl,'String',char(gui_state.prompt));
        set(answer_ctrl,'String','');
        set(button_ctrl,'String','Reveal');

    end

end


function [prompt_text,answer_text,context_text]=card_sides(card,mode)

    % Build the English side of the card
    eng_text=sprintf('%s [%s, %s]',char(card.english),...
                     char(card.source),char(card.form));

    % Choose the displayed side and the hidden answer
    if mode=="forward"
        prompt_text=string(eng_text);
        answer_text=card.hebrew;
        context_text="English -> Hebrew";
    else
        prompt_text=card.hebrew;
        answer_text=string(eng_text);
        context_text="Hebrew -> English";
    end

end


function font_name=select_font()

    % Define a cross-platform preference list with Hebrew coverage
    pref_fonts=["Noto Sans Hebrew","Arial Unicode MS","Arial",...
                "DejaVu Sans","Helvetica","SansSerif"];

    % Query installed fonts where the graphics system is available
    try
        avail_fonts=string(listfonts);
    catch
        avail_fonts=strings(0,1);
    end

    % Fall back to MATLAB's logical sans-serif font
    font_name="SansSerif";

    % Select the first preferred font installed on this system
    for n=1:numel(pref_fonts)
        if any(strcmpi(avail_fonts,pref_fonts(n)))
            font_name=pref_fonts(n);
            return
        end
    end

end


function grumble(mode,max_cards)
if ~isscalar(mode)||~ismember(mode,["forward","backward","both","gui"])
    error('mode must be ''forward'', ''backward'', ''both'', or ''gui''.');
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


