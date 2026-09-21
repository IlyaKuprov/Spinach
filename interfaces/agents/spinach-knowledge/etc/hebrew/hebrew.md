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
