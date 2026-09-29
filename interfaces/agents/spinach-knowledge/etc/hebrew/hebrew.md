# etc/hebrew/hebrew.m

## Call

```MATLAB
ncards = hebrew(mode, max_cards)
```

The source documents `hebrew()`, `hebrew(mode)`, `ncards=hebrew(mode,max_cards)`, and GUI mode; pass a string scalar such as `hebrew("gui")`. Defaults are `mode="both"` and `max_cards=Inf`. `mode` must be one of the exact lowercase string scalars `"forward"`, `"backward"`, `"both"`, or `"gui"`; the source compares these values directly. `max_cards` must be a non-negative numeric scalar. For example, `n = hebrew("both",0)` loads and counts the usable cards without starting a quiz.

## Workbook inputs and card construction

Place the vocabulary workbooks beside `hebrew.m`. The source names these files: `nouns.xlsx`, `adjectives.xlsx`, `verbs.xlsx`, `adverbs.xlsx`, `directions.xlsx`, `greetings.xlsx`, `languages.xlsx`, `numbers.xlsx`, `particles.xlsx`, `phrases.xlsx`, `prepositions.xlsx`, `pronouns.xlsx`, `proper_nouns.xlsx`, `quantifiers.xlsx`, `sex_terms.xlsx`, `weekdays.xlsx`, and `question_words.xlsx`. Each workbook's first column is the English prompt. Noun-like workbooks use English, masculine singular, feminine singular, masculine plural, feminine plural, invariant, and notes columns; adjectives omit invariant; verbs use English, infinitive, masculine singular, feminine singular, masculine plural, feminine plural, and notes.

The loader creates one card for each nonempty English/form pair only when that form contains at least one Hebrew Unicode character (U+0590–U+05FF). It records the source word class and grammatical-form name with each card. Empty input after filtering is an error. Thus preserve actual Hebrew Unicode in spreadsheet cells; this function does not translate, transliterate, or synthesise vocabulary. The MATLAB source itself contains no Hebrew vocabulary strings, so examples or translations cannot be recovered from it.

## Quiz behaviour and limits

- `"forward"` shows English first and reveals Hebrew after Enter; `"backward"` shows Hebrew first and reveals English. The English side includes source class and form.
- `"both"` chooses a direction randomly for each card. In console mode the cards are shuffled, then reshuffled after a full pass; an open-ended run ends with Ctrl+C. `max_cards` limits the console run.
- `"gui"` opens a one-button window: first click reveals the answer and the next advances. Direction is random, and the code avoids an immediate direction repeat when possible. The GUI is not capped by `max_cards`.

The workbook contents are runtime data external to this source file. The code filters for Hebrew-script code points but does not document or enforce a particular spreadsheet cell directionality/layout beyond preserving the text it reads.

Source: [hebrew.m in Spinach](https://github.com/IlyaKuprov/Spinach/blob/main/etc/hebrew/hebrew.m).
