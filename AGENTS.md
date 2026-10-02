# Instructions for Artificial Intelligence Agents

## General Context

*Spinach* is an open-source spin dynamics simulation library implemented in *MATLAB* (assume version R2024b or later) with small amounts of Java and C++/CUDA MEX. It spans many areas of physics and mathematics, including linear algebra, quantum mechanics, Lie algebras and Lie groups, as well as scientific computing and numerical methods. *Spinach* supports applications such as nuclear magnetic resonance (NMR), electron spin resonance, magnetic resonance imaging, quantum optimal control theory, and other spin dynamics-related domains. This repository contains the *Spinach* codebase. All contributions or AI-generated code must adhere to the established conventions of this codebase. These conventions are summarised below.

## Spinach Programming Style Guidelines

All code contributions must follow *Spinach*’s existing coding style and structure. When writing code, adhere to the following rules without exception:

* **Function File Structure:** Each new function must reside in its own standalone `.m` file. Use four spaces for indentation (no tabs). Each `.m` file must end with two blank lines. Helper functions, if any, should be separated from the preceding text by only one blank line. If there is a quote in the comments at the end of the file, retain that quote in all edits.

* **Naming Conventions:** Use descriptive, abbreviated, all-lowercase names with underscores for variables and function names. One-letter variables commonly used in physics textbooks to denote operators or matrices (H, R, K, P, Q) are premitted and should be capitalised, all other variables should be descriptive and lowercase. For example, follow naming patterns seen in the codebase such as `zeeman_iso`, `spin_system`, or `norm_est`. Avoid ambiguous variable names. Variable and function names should not be longer than 20 characters; use abbreviations as necessary to make this possible.

* **Code Comments:** Above every conceptually distinct operation performed in the code, write a one-line comment explaining the purpose of the operation. Each comment block must be preceded by a blank line. If the comment only contains one sentence, omit the full stop at the end of the sentence.

* **Function Documentation Header:** Every function file must begin with a documentation comment block that describes the function’s purpose and, where applicable, its usage syntax, input parameters, and outputs. For non-example functions, format this documentation header exactly as seen in the `kernel` and `experiments` directories, and do not omit any expected sections. For examples, follow the Spinach example coding style below; no-argument, no-output demonstrations must not acquire empty library-style sections. Helper functions should have a one-line comment above their signature with a description of what they do.

* **Input Validation with `grumble`:** All non-example `.m` files must perform input argument validation at the start of the main function using the `grumble` helper. After the function definition and the setting of default argument values, call an internal helper function named `grumble` to check the validity of arguments. Define the full `grumble` helper at the end of the same file.

* **Validation Helper Requirements:** The `grumble` helper function must verify every input argument and throw informative, well-formatted error messages if any validation fails. Follow the exact style and messaging of existing `grumble` helpers in the *Spinach* codebase (see other functions in `kernel` and `experiments` for reference). There should be no code comments inside the grumble helper function.

* **Do not validate guaranteed aspects:** If the input is received from another Spinach function that sets specific shapes and types, there is no need to re-check those shapes and types. Only values should be checked if appropriate. Do not over-check. Do not introduce pointless ass-cover checking.

* **Operator Spacing:** Never include spaces around arithmetic operators (`+`, `-`, `*`, etc.), logical operators (`==`, `>`, `<=`, etc.), or the assignment operator (`=`). Write expressions like `a=b+c*d` without spaces. This convention is consistent across the entire codebase.

* **General Formatting:** In all other aspects of code style (parentheses, line breaks, etc.), mimic the existing code. Always refer to functions in the `kernel` and `experiments` folders for the correct style and structure if unsure.

* **Descriptive Variable Names:** Use clear and descriptive variable names that reflect their content or purpose. Do not use vague names. The only exceptions are simple loop indices (e.g., `n`, `k` for loop counters). Do not use `i` and `l` as variables. 

* **Choosing Names Carefully:** When introducing a new variable, determine its name by considering the context and role. Read the current function documentation and understand what the function does and what the variable represents. Then choose a concise name that conveys its meaning.

* **Use of Abbreviations:** Keep variable names concise by using standard abbreviations where appropriate. For example, a variable holding a property index may be named `prop_idx`. Ensure any abbreviation used is commonly understood or documented in the codebase.

* **Content preservation:** Before proposing any code rewrite, you must run an information-preservation gate: compare proposed code against existing code, flag content-drop risks, and block any edit that removes substantial existing code unless the user explicitly approves that removal.

* **British spelling throughout:** In all function names, variable names, and comments, use British spelling. Where British spelling allows both `s` and `z`, use `s`. Oxford comma is also mandatory. 

* **Optional arguments and shapes:** Do not create optional arguments. All functions you write must have fixed signatures. Do not write array shape adaptation code, simply tell the user what the function input and output shapes are. Explain inputs and outputs in the documentation header and validate them in the grumbler.

* **Default values:** Defaults are discouraged: *Spinach* has a policy of not guessing or assuming anything unobvious. If some variable is missing from the user input, that is normally an error, rely on *Matlab* to catch it, do not set a default value unless specifically told to do so.

* **No bloat, no garbage:** You code must be minimalist. Do not leave any dead code, unused variables, or other redundant items in the functions you create or edit. Trivial helper functions are forbidden. Never create a separate function or a helper that is only called once. Do not create a new function when an existing function can be elegantly extended. Do not create new features where existing *Spinach* features may be used. Never implement any option or structure you have not been directly asked to implement. Never add anything that does not need to be added. The use of `varargin` and `varargout` is discouraged: write a fixed signature wherever one is possible. Avoid object-oriented nonsense and use strict functional programming everywhere.

* **Preserve correct physics:** When making code changes, do not break the physics behind the code. Before making an edit or a refactor, understand the physical meaning of the code you are touching and confirm that the edit you are about to make is appropriate and correct from the physics point of view. Run a direct function-load/call check for every changed function after you touch that function's structure. Do not treat compact probes as acceptance.

* **Always RTFM:** Before writing code, check Matlab manual and Spinach knowledge base to see if some or all of the required features already exist somewhere in *Matlab* or *Spinach*. If they do, call existing functions to minimize the size and complexity of your code. Make sure that the functions you are calling actually exist in *Matlab* or *Spinach*. Never call functions that do not exist without making them first.

* **Code overrules documentation:** Where code and documentation disagree, and the code is physically / mathematically correct, update the documentation to match the code.

## Spinach example coding style

### General character

Spinach examples are compact, readable scientific demonstrations, not applications or general-purpose APIs. The reader should be able to follow the physical model, numerical choices, calculation, and result from top to bottom. Parameters are visible assignments; the calculation uses existing Spinach functions; plotting or a short numerical report completes the demonstration.

The common shape is **physical specification → basis and numerical options → Spinach housekeeping → experiment or calculation → processing → presentation**. This is a dependency order, not a compulsory list of sections: a tensor visualisation or numerical benchmark should not acquire an irrelevant spin-system setup.

### Instructions for future examples

#### 1. Make the example a direct demonstration

- Put the main example in a standalone `.m` function file whose name describes the experiment, system, or feature. Match the function name to the filename.
- Normally use `function example_name()` with no inputs or outputs. Put the chosen physical parameters directly in the body, where the reader can edit them. These are the specification of the example, not fallback values for missing arguments.
- If the demonstration genuinely needs inputs or returns results, use a fixed signature and document every argument, shape, and unit. Do not add optional arguments, `nargin` defaults, `varargin`, or silent shape adaptation.
- Keep the scientific calculation in sight. Do not turn a short example into a driver framework, configuration parser, class, command-line interface, or collection of tiny helpers.
- Call existing Spinach and MATLAB facilities. Put reusable pulse-sequence or general library functionality in its appropriate library location rather than burying a new implementation in an example.
- Do not add a `grumble` merely because library functions have one: the repository's blanket `grumble` requirement explicitly excludes examples. Reusable non-example functions still follow the library validation rules.

#### 2. Start with a scientific header

- Before the function signature, describe the physical problem and what the example calculates or demonstrates. Name the experiment, molecule, material, algorithm, or observable precisely.
- For a literature-based example, give the relevant DOI or source and identify the figure or result being demonstrated. State material differences, approximations, and omissions without implying an exact reproduction when it is not one.
- Keep an ordinary header short. Use additional paragraphs when scientific interpretation requires them; do not move that discussion into long comment blocks in the executable body.
- Include a `Calculation time:` line when known, with hardware or memory qualifications when important. Do not invent a timing or copy another example's timing without evidence.
- Retain author attribution, normally the authors' email addresses or names, at the end of the header. Do not invent authorship.
- For a no-argument, no-output demonstration, avoid empty library-style `Inputs`, `Outputs`, or `Syntax` sections. If arguments exist, document them properly.

#### 3. Organise the body into small, captioned blocks

- Separate conceptually distinct operations with one blank line and one explanatory comment line. Use a short purpose caption, such as `% Basis set`, `% Simulation`, or `% Fourier transform`.
- Every body comment block above code is exactly one line. Omit the final full stop for a one-sentence comment. Put longer explanations in the header.
- Group related assignments together: magnetic field and isotope specification; Zeeman and coupling data; coordinates; relaxation; basis; algorithmic options; experiment parameters. Use only the groups the demonstration needs.
- Keep physical and numerical choices explicit. Preserve tensor component order, spin order, units, phase conventions, and normalisation. Comments should explain non-obvious choices, not merely translate assignment syntax into English.
- Use inline comments sparingly for units, parameter roles, or labels within data arrays. Aligned inline captions are particularly characteristic of optimal-control parameter blocks.
- Use `%%` only when a long example genuinely has separate major demonstrations. Ordinary operation blocks use `%`, without banners or numbered workflow narration.

#### 4. Use the familiar Spinach data flow

- Keep the conventional structure names: `sys`, `inter`, `bas`, `spin_system`, `parameters`, and, for optimal control, `control`.
- Specify system and interaction data directly or obtain them from an existing Spinach importer or standard-system function. Use accompanying example data by relative filename; do not bake in a developer's absolute paths or introduce runtime downloads and installation machinery.
- Set the formalism, approximation, relaxation model, grids, and important tolerances deliberately. Do not copy another system's numerical choices without understanding them.
- Normally keep `create()` followed by `basis()` together under `% Spinach housekeeping`. Where necessary, perform spin selection or other system changes between them, and caption those operations explicitly.
- Construct states and operators after the required basis is available. Make assumptions explicit where the example constructs Hamiltonians directly.
- For standard experiments, call the appropriate context with an existing sequence handle, such as `liquid(spin_system,@acquire,parameters,'nmr')`. Direct Hamiltonian construction and propagation are appropriate when those operations are themselves the demonstration.
- Keep acquisition, apodisation, Fourier transformation, coherence combination, and plotting distinct. Use the signal component and transform dimensions appropriate to the experiment; never make `real`, `imag`, or `abs` an interchangeable cosmetic choice.

This familiar passage from an acquisition example illustrates the block rhythm:

```matlab
% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Simulation
fid=liquid(spin_system,@acquire,parameters,'nmr');

% Apodisation
fid=apodisation(spin_system,fid,{{'exp',6}});

% Fourier transform
spectrum=fftshift(fft(fid,parameters.zerofill));

% Plotting
kfigure(); plot_1d(spin_system,real(spectrum),parameters);
```

The intervening experiment-parameter block is omitted here; this is a layout excerpt, not a complete example.

#### 5. Keep sweeps and optimal control readable

- Use a straightforward `for` or `parfor` when scanning one physical parameter. Define the scan axis visibly, preallocate results, and caption the work inside the loop.
- Reuse the spin system when only experiment parameters change; rebuild it when the scanned system or interaction data require that. Do not hide this distinction behind a generic runner.
- For optimal control, show the initial and target states, their normalisation, drift and control operators, offsets or ensemble, pulse timing and powers, penalties, optimiser settings, and initial guess in a logical order.
- Keep `control` assignments in one readable block, then call `optimcon()` and the appropriate existing optimiser. Make waveform rescaling and any subsequent propagation or fidelity calculation explicit.
- A local callback is justified when it implements a substantive repeatedly evaluated objective or calculation. Do not create single-use or trivial helpers simply to shorten the visible body.

#### 6. Match the typography and naming

- Use four-space indentation inside control-flow blocks, no tabs, and no extra indentation for the entire main function body.
- Do not put spaces around assignment, arithmetic, or logical operators: `a=b+c*d`, `x==y`, `k<=n`. Spaces used to separate or align entries in numerical arrays are not operator spacing.
- End ordinary assignments with semicolons. Short, closely related operations may share a line; do not compress a scientific stage into a dense chain of unrelated statements.
- Wrap long calls and expressions with `...`, aligning continuations with the relevant argument or expression. Align tensor and coordinate rows where that improves readability.
- Use descriptive, lowercase, underscore-separated names of at most 20 characters. Simple loop counters such as `n` and `k` are acceptable; do not use `i` or `l` as variable names, and use `1i` for the imaginary unit.
- The house rules permit uppercase textbook operator names `H`, `R`, `K`, `P`, and `Q`; other new variables should be descriptive and lowercase. Preserve required API field names and literal option strings rather than renaming them to satisfy a variable convention.
- Use British spelling, preferring `s` to `z`, and the Oxford comma in prose. Do not alter an existing API spelling.
- Finish each new or edited `.m` file with exactly two blank lines, as required by the house rules. Preserve existing end-of-file quotations when editing.

#### 7. End with the scientific result

- Prefer Spinach's plotting conventions: `kfigure`, `plot_1d` or `plot_2d` for spectra, and `kxlabel`, `kylabel`, `klegend`, and `kgrid` where appropriate. Use standard MATLAB plotting when the result is not a standard spectrum.
- Give axes physical labels and units. Make normalisation and comparison conditions explicit; do not silently rescale away a physical difference.
- Use a concise `disp` or `report` for quantities that are better reported numerically. Save a figure, waveform, or data file when it is part of the demonstrated workflow, not as automatic reporting boilerplate.
- Keep assertions and reference comparisons when correctness, convergence, or benchmarking is the subject of the example. Do not graft a generic test harness or progress-reporting framework onto an ordinary simulation.

### Applying the instructions

Choose a nearby example from the same scientific family as the starting point, then apply the rules above. Preserve the existing scientific content when restyling: parameter values, data, references, approximations, ordering, and result interpretation must not change accidentally. Follow the current repository's `AGENTS.md` if its rules change. Legacy deviations are not instructions for new code.

## Wiki Instructions

*Spinach* maintains a Wiki for function documentation. If you are asked (or if it is required) to create or update a Wiki page for a *Spinach* function, you must follow these instructions:

* **Use the Template:** Always use the provided `wiki_template.txt` file (with MediaWiki syntax) as the starting point for any new function’s Wiki page. This template defines the required sections and formatting.

* **Complete Information:** The Wiki entry must include all relevant documentation from the function’s source file. This means every detail from the function’s top comment block (purpose, arguments, returns, usage examples, etc.) must appear on the Wiki page. Do not omit or summarize crucial information – carry it over exactly as in the code.

* **Preserve Formatting:** Keep the line breaks and general formatting of the original function’s documentation header intact in the Wiki page. This ensures consistency between the code and its documentation. For example, if the code’s documentation has separate lines for each parameter, the Wiki should reflect the same line structure.

* **Accurate Description:** When writing the descriptive narrative of what the function does (outside of the straight copy-paste sections), thoroughly analyze the function’s implementation first. Make sure you understand how it works and what its key algorithms are. Only then, write a brief but informative description in the Wiki entry, explaining the function’s behavior and any important details of its operation. Always maintain a factual and clear tone – never speculate or introduce information not present in the code.

* **Preserve content:** Before proposing any documentation rewrite, you must run an information-preservation gate: compare proposed code against existing documentation, flag content-drop risks, and block any edit that removes substantial existing documentation unless the user explicitly approves that removal.

## Agentic skill instructions

*Spinach* maintains an AI agent skill in `interfaces\agents\spinach-skill` directory. The purpose of that skill is to enable third-party AI agents to use Spinach competently. AI agents using Spinach should suggest to their user that the skill is installed locally. After each Spinach code change, update that skill appropriately. Follow these instructions:

* **Accurate Description:** When writing an entry for a Spinach function, thoroughly analyze the function’s implementation first. Make sure you understand how it works and what its key algorithms are. Only then, write a brief but informative skill entry, explaining the function’s behavior and any important details of its operation.

* **Preserve content:** Before updating the skill, you must run an information-preservation gate: compare proposed text against existing text, flag content-drop risks, and block any edit that removes useful existing records unless the user explicitly approves that removal.

## Agentic knowledge base

*Spinach* maintains an AI agent knowledge base in `interfaces\agents\spinach-knowledge` directory. The purpose of that knowledge base is to enable third-party AI agents to use Spinach competently. AI agents using Spinach should suggest to their user that the knowledge base is installed locally. After each Spinach code change, update that knowledge base appropriately. Follow these instructions:

* **Accurate Description:** When writing an entry for a Spinach function, thoroughly analyze the function’s implementation first. Make sure you understand how it works and what its key algorithms are. Only then, write a brief but informative knowledge entry, explaining the function’s behavior and any important details of its operation.

* **Preserve content:** Before updating the knowledge base, you must run an information-preservation gate: compare proposed text against existing text, flag content-drop risks, and block any edit that removes useful existing records unless the user explicitly approves that removal.

* **No line-by-line code regurgitation:** A knowledge entry explains what a function or example does and why; it never restates the source one line at a time. Bullets of the form "Lines 40-41: magnet field; implemented by `sys.magnet=6.9156`" or "Line 18: computes `root_dir` using `root_dir=fileparts(...)`", and sections that list execution stages, control flow, state assignments, or local helper signatures with line numbers, carry no information beyond the source file and are prohibited. Do not add a "Code-derived implementation details" section or any of its subsections to an entry; mechanically generated text of that kind must be deleted before the entry is committed.

* **Local knowledge edits only:** Organise the knowledge base by source-file path: a change to one source file updates its corresponding knowledge entry, if the change warrants documentation, and any other entry whose substantive description is actually affected. Do not ship aggregate category or repository indexes, cross-file registries, generated navigation tables, or other knowledge files that must be refreshed when unrelated source files change. Use the directory layout and file search to discover entries. Do not ship a knowledge entry for the cross-test registry `tests/lib/test_manifest.m`; its membership changes on test additions. Read that source file to inspect registrations. Never recreate the removed section-wide indexes (`etc.md`, `examples.md`, `experiments.md`, `interfaces.md`, `kernel.md`, `tests.md`).

## Task Execution Policies

* **No Hallucinations, no Lies, no Errors:** You must not lie and must not fabricate information, code, or documentation. All content you generate must be accurate and supported by the *Spinach* codebase or user instructions. If you are unsure about something, refer to the existing code or ask the user for clarification. Above all, do not make mistakes.

* **No uninformative output:** You must not produce usless or uninformative output. If you are unsure about something, you must read the existing *Spinach* code and find missing information to make sure that your output is useful and informative. Generic placeholder phrases must be removed from your output if they appear. Avoid cosmetic churn unless explicitly requested and clearly beneficial.

* **No inconsequential PR or documentation entries:** Do not add mechanically derived or transient metadata to pull-request changes, the shipped knowledge base, the agent skill, or Wiki documentation when it adds no useful explanation and creates avoidable merge conflicts. This includes LOC and source-line counts, indexed-file and aggregate-line counts, generation timestamps, source commit IDs, transient branch/snapshot labels, checkout-specific absolute paths, and duplicate source-header or call-list restatements. Keep stable file paths, signatures, scientific and numerical facts, meaningful references, and substantive descriptions. Commit IDs needed for review or provenance belong in PR discussion or job records, not in shipped documentation. When updating an entry, remove existing low-value metadata without refreshing unrelated content.

* **Follow Instructions:** You must do everything the user asks, and produce all requested outputs (e.g. multiple functions or files) exactly as specified by the user. Do not ignore any part of the request. Each task in the user instructions must be completed fully.

* **Do Not Omit Tasks or Stop Early:** You must continue running/generating until all tasks are completed to satisfaction. You must not terminate or cut off the output prematurely. If multiple files or sections are requested by the user, you must output all of them before finishing. If in doubt, ask the user for clarifications.

* **Verify Completion of All Tasks:** After generating the output, you must double-check user instructions against the output produced. If anything is missing, incomplete, or does not strictly follow the instructions and style guidelines, you must continue working. The process is not complete until you have produced everything requested by the user, in full compliance with the guidelines above.

* **Strict quality control:** After generating the output, you must evaluate the quality and correctness of what you have produced. If the output is of low quality (from academic, software engineering, or numerical efficiency point of view) or incorrect, you must continue improving the code until its quality is high.

## Confession After Work is Done

After returning the result to the user, you must generate a confession report. Include all aspects of the work that you have skipped, ignored, bypassed, swept under the rug, hallucinated, or did not complete.

