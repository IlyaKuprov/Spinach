# etc/hebrew/hebrew.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/etc/hebrew/hebrew.m`
- Signature: `ncards=hebrew(mode,max_cards) % #NWIKI #NHEAD`
- Total lines: 389

## Purpose

IK's Hebrew flashcards function. The Excel files should contain Hebrew vocabulary in separate spreadsheets: nouns.xlsx -English, masculine singular, feminine singular, masculine plural, feminine plural, invariant, notes adverbs.xlsx, directions.xlsx, greetings.xlsx, languages.xlsx, numbers.xlsx, particles.xlsx, phrases.xlsx, prepositions.xlsx, pronouns.xlsx, proper_nouns.xlsx, quantifiers.xlsx, sex_terms.xlsx, weekda

## Physical / mathematical content

- This file belongs to the `etc` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `load_cards()`, `read_cards()`, `clean_text()`, `contains_hebrew()`, `sprintf()`, `select_font()`, `button_press()`, `randi()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 42-43: Set default direction mode; implemented by `if nargin<1`.
- Lines 47-48: Set default number of cards; implemented by `if nargin<2`.
- Lines 52-53: Normalise and validate inputs; implemented by `mode=lower(string(mode))`.
- Lines 56-57: Locate the spreadsheet directory; implemented by `this_file=mfilename('fullpath')`.
- Lines 60-61: Load all usable flashcards; implemented by `cards=load_cards(root_dir)`.
- Lines 64-65: Return immediately for spreadsheet-loading tests; implemented by `if max_cards==0`.
- Lines 69-70: Start the graphical flashcard interface; implemented by `if mode=="gui"`.
- Lines 75-76: Start with a random card order; implemented by `card_order=randperm(ncards)`.
- Lines 79-80: Main flash card loop; implemented by `while shown<max_cards`.
- Lines 82-83: Reshuffle after all cards have been used; implemented by `if order_pos>ncards`.
- Lines 89-90: Select the next card; implemented by `card_idx=card_order(order_pos)`.
- Lines 94-95: Randomise direction in mixed mode; implemented by `card_mode=mode`.
- Lines 104-105: Show the selected card; implemented by `show_card(cards(card_idx,:),card_mode)`.

### Control flow inferred from the code

- Line 43: conditional branch on `nargin<1`.
- Line 48: conditional branch on `nargin<2`.
- Line 65: conditional branch on `max_cards==0`.
- Line 70: conditional branch on `mode=="gui"`.
- Line 80: `while` loop over `shown<max_cards`.
- Line 83: conditional branch on `order_pos>ncards`.
- Line 96: conditional branch on `mode=="both"`.
- Line 97: conditional branch on `rand()<0.5`.

### Key state/data transformations

- Lines 44: computes `mode` using `mode="both"`.
- Lines 49: computes `max_cards` using `max_cards=inf`.
- Lines 57: computes `this_file` using `this_file=mfilename('fullpath')`.
- Lines 58: computes `root_dir` using `root_dir=fileparts(this_file)`.
- Lines 61: computes `cards` using `cards=load_cards(root_dir)`.
- Lines 62: computes `ncards` using `ncards=height(cards)`.
- Lines 76: computes `card_order` using `card_order=randperm(ncards)`.
- Lines 77: computes `order_pos` using `order_pos=1; shown=0`.
- Lines 90: computes `card_idx` using `card_idx=card_order(order_pos)`.
- Lines 92: computes `shown` using `shown=shown+1`.
- Lines 95: computes `card_mode` using `card_mode=mode`.

### Local helper functions

- Line 111: `load_cards()` — `function cards=load_cards(root_dir)`. Start from an empty table
  - Representative operation: `cards=table(strings(0,1),strings(0,1),strings(0,1),strings(0,1), 'VariableNames',{'english','hebrew','form','source'})`.
  - Representative operation: `'VariableNames',{'english','hebrew','form','source'})`.
- Line 164: `read_cards()` — `function cards=read_cards(file_name,source,form_names,form_cols)`. Read the spreadsheet as strings
  - Representative operation: `T=readtable(file_name,'TextType','string','PreserveVariableNames',true)`.
  - Representative operation: `english=clean_text(T{:,1})`.
- Line 203: `clean_text()` — `function text=clean_text(text)`. Convert spreadsheet cells into trimmed strings
  - Representative operation: `text=string(text)`.
  - Representative operation: `text(ismissing(text))=""`.
- Line 212: `contains_hebrew()` — `function tf=contains_hebrew(text)`. Detect Hebrew code points
  - Representative operation: `tf=~isempty(regexp(char(text),'[\x{0590}-\x{05FF}]','once'))`.
- Line 219: `show_card()` — `function show_card(card,mode)`. Build the English side of the card
  - Representative operation: `eng_text=sprintf('%s [%s, %s]',char(card.english), char(card.source),char(card.form))`.
  - Representative operation: `char(card.source),char(card.form))`.
- Line 239: `show_gui()` — `function show_gui(cards)`. Select a Hebrew-capable display font
  - Representative operation: `gui_font=select_font()`.
  - Representative operation: `gui_fig=figure('Name','Hebrew flashcards','NumberTitle','off', 'MenuBar','none','ToolBar','none','Color',[1 1 1], 'Units','pixels','Position',[100 100 760 430])`.
- Line 285: `button_press()` — `function button_press(~,~)`. Reveal or advance depending on the current state
  - Representative operation: `if gui_state.revealed`.
  - Representative operation: `next_card()`.
- Line 298: `next_card()` — `function next_card()`. Select a random card without immediate repetition where possible
  - Representative operation: `new_idx=randi(height(gui_state.cards))`.
  - Representative operation: `if height(gui_state.cards)>1&&new_idx==gui_state.card_idx`.

## Syntax

```matlab
hebrew()
hebrew(mode)
ncards=hebrew(mode,max_cards)
hebrew('gui')
The mode may be 'forward', 'backward', 'both', or 'gui'. In forward mode,
English is shown first and Hebrew is revealed after <Enter>. In
backward mode, Hebrew is shown first and English is revealed after
<Enter>. In both mode, the direction is random on every card. Exit
an open-ended run with CTRL+C. In gui mode, a one-button graphical
flashcard window is opened; the button reveals the answer on the
first click and advances to another randomly selected card on the
second click.
```

## Implementation structure

- IK's Hebrew flashcards function. The Excel files should contain
- Hebrew vocabulary in separate spreadsheets:
- nouns.xlsx -English, masculine singular, feminine singular,
- masculine plural, feminine plural, invariant,
- notes
- adverbs.xlsx, directions.xlsx, greetings.xlsx, languages.xlsx,
- numbers.xlsx, particles.xlsx, phrases.xlsx,
- prepositions.xlsx, pronouns.xlsx,
- proper_nouns.xlsx, quantifiers.xlsx,
- sex_terms.xlsx, weekdays.xlsx -same columns as
- nouns.xlsx
- adjectives.xlsx -English, masculine singular, feminine singular,

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `lower()`, `string()`, `grumble()`, `mfilename()`, `fileparts()`, `load_cards()`, `height()`, `show_gui()`, `randperm()`, `card_order()`, `show_card()`, `cards()`, `table()`, `strings()`, `fullfile()`, `read_cards()`.
