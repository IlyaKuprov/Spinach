# etc/hebrew/hebrew.m

- Signature: `ncards=hebrew(mode,max_cards) % #NWIKI #NHEAD`

## Purpose

IK's Hebrew flashcards function. The Excel files should contain Hebrew vocabulary in separate spreadsheets: nouns.xlsx -English, masculine singular, feminine singular, masculine plural, feminine plural, invariant, notes adverbs.xlsx, directions.xlsx, greetings.xlsx, languages.xlsx, numbers.xlsx, particles.xlsx, phrases.xlsx, prepositions.xlsx, pronouns.xlsx, proper_nouns.xlsx, quantifiers.xlsx, sex_terms.xlsx, weekda

## Physical / mathematical content

## Numerical / algorithmic content

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
