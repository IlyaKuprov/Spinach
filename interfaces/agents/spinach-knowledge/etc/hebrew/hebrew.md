# etc/hebrew/hebrew.m

- Signature: `ncards=hebrew(mode,max_cards) % #NWIKI #NHEAD`

Loads Hebrew vocabulary flashcards from Excel spreadsheets. Defaults are `mode="both"` and `max_cards=Inf`; mode names are case-insensitive. `max_cards=0` loads the spreadsheets and returns the number of cards without presenting them.

## Spreadsheet format

These files use columns for English, masculine singular, feminine singular, masculine plural, feminine plural, invariant, and notes: `nouns.xlsx`, `adverbs.xlsx`, `directions.xlsx`, `greetings.xlsx`, `languages.xlsx`, `numbers.xlsx`, `particles.xlsx`, `phrases.xlsx`, `prepositions.xlsx`, `pronouns.xlsx`, `proper_nouns.xlsx`, `quantifiers.xlsx`, `sex_terms.xlsx`, `weekdays.xlsx`, and `question_words.xlsx`.

`adjectives.xlsx` uses English, the four gender/number forms, and notes. `verbs.xlsx` uses English, infinitive, the four gender/number forms, and notes. Nonempty English/Hebrew forms containing Hebrew characters become cards; notes are not cards. An empty card set is an error.

## Modes

- `forward`: show English, then reveal Hebrew after Enter.
- `backward`: show Hebrew first, then reveal English.
- `both`: choose a direction randomly for each card.
- `gui`: open a one-button window; click to reveal the answer, then click again to advance. Direction is random, with immediate repeats avoided when possible.

Console cards are shuffled and reshuffled after all have been shown. `max_cards` limits the console run. Press Ctrl+C to stop an open-ended run.