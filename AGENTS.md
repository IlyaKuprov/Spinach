# Instructions for Artificial Intelligence Agents

## General Context

*Spinach* is an open-source spin dynamics simulation library implemented in *MATLAB* (assume version R2026a or later) with small amounts of Java and C++/CUDA MEX. It spans many areas of physics and mathematics, including linear algebra, quantum mechanics, Lie algebras and Lie groups, as well as scientific computing and numerical methods. *Spinach* supports applications such as nuclear magnetic resonance, electron spin resonance, magnetic resonance imaging, quantum optimal control theory, and other spin dynamics related domains. This repository contains the *Spinach* codebase. All contributions or AI-generated code must adhere to the established conventions of this codebase. These conventions are summarised below.

## Task Execution Policies

### Scope, Accuracy, and Communication

* **Follow Instructions and Verify Completion:** Do everything the user asks, producing every requested file, function, section, and other output exactly as specified. Do not ignore any part of the request, cut off the output, or stop early. If in doubt, ask the user for clarification. Before finishing, double-check the output against every user instruction and the guidelines above, and evaluate its correctness and academic, software-engineering, and numerical-efficiency quality. Continue working until nothing is missing, incomplete, incorrect, non-compliant, or low-quality.

* **No Hallucinations, no Lies, no Errors:** You must not lie and must not fabricate information, code, or documentation. All content you generate must be accurate and supported by the *Spinach* codebase or user instructions. If you are unsure about something, refer to the existing code or ask the user for clarification. Above all, do not make mistakes.

* **No uninformative output:** You must not produce usless or uninformative output. If you are unsure about something, you must read the existing *Spinach* code and find missing information to make sure that your output is useful and informative. Generic placeholder phrases must be removed from your output if they appear. Avoid cosmetic churn unless explicitly requested and clearly beneficial.

* **Pull request merging requires permission:** When asked to submit a pull request, prepare it and address Codex bot review comments, but do not merge the pull request. Merging is something the user will do after reviewing the pull request.

### Preservation and Scientific Correctness

* **Preserve content:** Before proposing a code or documentation rewrite or updating the agent skill or knowledge base, run an information-preservation gate: compare proposed content against existing content and flag content-drop risks. Block removal of substantial existing code or documentation, or useful existing skill or knowledge-base records, unless the user explicitly approves that removal.

* **Preserve correct physics:** When making code changes, do not break the physics behind the code. Before making an edit or a refactor, understand the physical meaning of the code you are touching and confirm that the edit you are about to make is appropriate and correct from the physics point of view. Run a direct function-load/call check for every changed function after you touch that function's structure. Do not treat compact probes as acceptance.

* **Code overrules documentation:** Where code and documentation disagree, and the code is physically / mathematically correct, update the documentation to match the code.

## Spinach Programming Style Guidelines

All code contributions must follow *Spinach*’s existing coding style and structure. When writing code, adhere to the following rules:

### Planning and Reuse

* **Always RTFM:** Before writing code, check Matlab manual and Spinach knowledge base to see if some or all of the required features already exist somewhere in *Matlab* or *Spinach*. If they do, call existing functions to minimize the size and complexity of your code. Make sure that the functions you are calling actually exist in *Matlab* or *Spinach*. Never call functions that do not exist without making them first.

* **No bloat, no garbage:** Your code must be minimalist. Remove dead code, unused variables, redundant conversions, unnecessary aliases, speculative branches, unrequested options, and other redundant items from functions you create or edit. Do not create trivial or single-use helpers; the required local `grumble` and substantive repeatedly evaluated numerical callbacks, objectives, or calculations are explicit exceptions. Elegantly extend an existing function or reuse existing Spinach features instead of adding a parallel implementation. Never implement an option or structure you have not been directly asked to implement, and never add anything unnecessary. Avoid object-oriented nonsense and use strict functional programming everywhere.

### Interfaces and Defaults

* **Optional arguments and shapes:** Prefer fixed signatures; avoid optional arguments, `varargin`, and `varargout` wherever a fixed signature is possible. For new inputs, require all arguments: do not add `nargin` defaults, optional arguments, `varargin`, or `varargout`. Do not write array shape adaptation code. Document input and output shapes, units, and meanings in the header and validate them in the grumbler where required. Preserve existing documented call patterns when modifying a legacy function; do not introduce unsolicited API breaks.

* **Default values:** Defaults are discouraged: *Spinach* has a policy of not guessing or assuming anything unobvious. If some variable is missing from the user input, that is normally an error, rely on *Matlab* to catch it, do not set a default value unless specifically told to do so.

### Input Validation

* **Input Validation with `grumble`:** All non-example `.m` files must perform input argument validation at the start of the main function using the `grumble` helper. After the function definition and the setting of default argument values, call an internal helper function named `grumble` to check the validity of arguments. Define the full `grumble` helper at the end of the same file. Do not add a `grumble` to an example merely because library functions have one; reusable non-example functions still require it.

* **Validation Helper Requirements:** Put ordinary input checks in the local `grumble`, verifying every input argument’s caller-controlled properties needed by the algorithm. Follow the exact style and messaging of existing helpers in `kernel` and `experiments`: concise, informative, well-formatted errors, and no comments inside the helper. Keep genuinely computed-domain checks beside the quantity they inspect, as in relaxation calibration. Reject unsupported cases clearly; do not catch an error and silently substitute different physics.

* **Do not validate guaranteed aspects:** Do not recursively revalidate a Spinach structure or recheck properties already guaranteed by its producer. If another Spinach function sets specific shapes and types, check only input values where appropriate. Do not over-check or introduce pointless ass-cover checking.

### File Layout, Naming, and Formatting

* **Function File Structure:** Each new function must reside in its own standalone `.m` file. Use four spaces for indentation (no tabs). Each `.m` file must end with exactly two blank lines. Helper functions, if any, should be separated from the preceding text by only one blank line. If there is a quote in the comments at the end of the file, retain that quote in all edits. Every `end` keyword must be correctly indented. Maintain formatting symmetry between opening and closing keywords: when its opening statement is followed by a blank line, the closing `end` must be preceded by a blank line. If a Matlab command or an array is broken up into multiple lines using `...`, the lines must be logically and elegantly indented: decimal point alignment for numerical arrays, structural alignment for multiple arguments, etc.

* **Naming Conventions:** Use descriptive, concise, all-lowercase, underscore-separated variable and function names of at most 20 characters. Read the current function documentation and consider each variable’s context, content, and role before naming it; avoid ambiguous or vague names. Use standard, commonly understood or documented abbreviations, such as `prop_idx` for a property index. Follow nearby naming patterns when in doubt. Uppercase textbook operator names (`H`, `R`, `K`, `P`, `Q`) and matrix names matching their mathematical notation are permitted. Loop counters should be single lowercase letters, such as `n` and `k`; do not use `i` or `l` as variables, and use `1i` for the imaginary unit. Use British spelling in function names, variable names, and comments, preferring `s` to `z` where both are allowed, and use the Oxford comma. Preserve required API field names, literal option strings, and existing API spellings rather than renaming them to satisfy these conventions.

* **Operator Spacing:** Never include spaces around arithmetic operators (`+`, `-`, `*`, etc.), logical operators (`==`, `>`, `<=`, etc.), or the assignment operator (`=`). Write expressions like `a=b+c*d` without spaces. Spaces used to separate or align entries in numerical arrays are not operator spacing.

* **General Formatting:** When in doubt about formatting, mimic the existing code. Refer to functions in the `kernel` and `experiments` folders for the correct style and structure if unsure.

### Comments and Function Headers

* **Code Comments:** Above every conceptually distinct operation, write exactly one line explaining its purpose, preceded by one blank line. Omit the final full stop for a one-sentence comment. Put extended explanations in the function header rather than long body comment blocks.

* **Function Documentation Header:** Every function file must begin with a documentation comment block that describes the function’s purpose and, where applicable, its usage syntax, input parameters, and outputs. For non-example functions, format this documentation header exactly as seen in the `kernel` and `experiments` directories, and do not omit any expected sections. For examples, follow the Spinach example coding style below; no-argument, no-output demonstrations must not acquire empty library-style `Inputs`, `Outputs`, or `Syntax` sections. Helper functions should have a one-line comment above their signature with a description of what they do. Preserve existing author attribution.

## Spinach Kernel Coding Elegance

Keep the scientific model, algorithm, and numerical costs directly readable. Apply these rules alongside the programming style above; retain established API behaviour. Read this section before writing or refactoring a Spinach function, and apply the final review before delivery.

### Write the calculation, not a framework

1. **Start from a nearby kernel analogue.** Read its scientific contract and implementation before designing the change; apply the existing-capability and minimalism rules above.

2. **Keep the mathematics visible.** Use expressions that a domain scientist can match to the defining equation. Keep small fixed transforms explicit when that is clearer than a generic implementation. Introduce an intermediate when it names a physical quantity, exposes a numerical stage, or avoids meaningful repeated work—not merely to rename an expression.

3. **Use direct data flow.** Arrange the body in dependency order, with one conceptually distinct operation per captioned block. Prefer ordinary matrices, arrays, cells, and existing Spinach structures. Do not introduce a class, dispatcher framework, configuration parser, or mutable cache for a routine calculation.

4. **Abstract scientific reuse, not typing effort.** Apply the helper and reuse rules above. Do not force distinct physical cases through a more complicated universal mechanism merely to eliminate repeated lines.

### Make the numerical representation earn its cost

5. **Compute only what the caller needs.** If an operator action suffices, do not materialise the full operator or exponential. Do not allocate a trajectory, Hessian, or accumulated propagator unless requested. Preserve sparse, factorised, or matrix-free structure until an operation genuinely requires expansion.

6. **Use MATLAB operations where they express the mathematics clearly.** Prefer direct array algebra and established sparse constructors to hand-built book-keeping. Keep loops when they naturally traverse spins, tensor ranks, pulse slices, or independent cases. Do not replace a readable loop with dense broadcasting, elaborate indexing, or a large temporary merely to call the result vectorised. Preallocate growing results and avoid repeating expensive conversions inside loops.

7. **Give numerical decisions a reason.** Derive scales, convergence decisions, and bounds from the algorithm and data where possible; reuse the relevant Spinach tolerance policy. Distinguish exact mathematical constants from numerical tolerances and performance cutoffs. Do not introduce an unexplained epsilon, iteration cap, size threshold, or regulariser. Do not silently change established thresholds during restyling, and do not claim a performance improvement without measurement.

8. **Preserve the scientific contract.** Keep units, signs, normalisation, basis order, shapes, and operator conventions explicit. For propagation, retain the convention `exp(-1i*L*t)`. Do not silently symmetrise, renormalise, clip, transpose, or regularise an input to make the calculation work; such operations need a documented mathematical role and compatibility with the requested behaviour.

### Review before delivery

9. **Review for subtraction before delivery.** Apply the minimalism rules above while keeping the scientific calculation locally understandable; compactness is not code golf. Compare against the original to preserve scientific content and established behaviour. For every changed MATLAB file, run the required house-style checker, MATLAB `checkcode`, and the validation required by the repository and change; style compliance does not prove numerical correctness.

**Final review question:** Can a Spinach scientist see the model, the algorithm, and the numerical costs directly, without mentally dismantling software scaffolding? If not, simplify the design before polishing its formatting.

## Spinach example coding style

### General character

Spinach examples are compact, readable scientific demonstrations, not applications or general-purpose APIs. The reader should be able to follow the physical model, numerical choices, calculation, and result from top to bottom. Parameters are visible assignments; the calculation uses existing Spinach functions; plotting or a short numerical report completes the demonstration.

The common shape is **physical specification → basis and numerical options → Spinach housekeeping → experiment or calculation → processing → presentation**. This is a dependency order, not a compulsory list of sections: a tensor visualisation or numerical benchmark should not acquire an irrelevant spin-system setup.

### Instructions for future examples

#### 1. Make the example a direct demonstration

- Put the main example in a standalone `.m` function file whose name describes the experiment, system, or feature. Match the function name to the filename.

- Normally use `function example_name()` with no inputs or outputs. Put the chosen physical parameters directly in the body, where the reader can edit them. These are the specification of the example, not fallback values for missing arguments.

- If the demonstration genuinely needs inputs or returns results, use a fixed signature and document every argument, shape, and unit, following the interface rules above.

- Keep the scientific calculation in sight. Do not turn a short example into a driver framework, configuration parser, class, command-line interface, or collection of tiny helpers.

- Put reusable pulse-sequence or general library functionality in its appropriate library location rather than burying a new implementation in an example.

#### 2. Start with a scientific header

- Before the function signature, describe the physical problem and what the example calculates or demonstrates. Name the experiment, molecule, material, algorithm, or observable precisely.

- For a literature-based example, give the relevant DOI or source and identify the figure or result being demonstrated. State material differences, approximations, and omissions without implying an exact reproduction when it is not one.

- Keep an ordinary header short. Use additional paragraphs when scientific interpretation requires them; do not move that discussion into long comment blocks in the executable body.

- Include a `Calculation time:` line when known, with hardware or memory qualifications when important. Do not invent a timing or copy another example's timing without evidence.

- Retain author attribution, normally the authors' email addresses or names, at the end of the header. Do not invent authorship.

#### 3. Organise the body into small, captioned blocks

- Follow the one-line body-comment rule above, using short purpose captions such as `% Basis set`, `% Simulation`, or `% Fourier transform`.

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

#### 6. Match the typography and naming

- Apply the programming style above, with no extra indentation for the entire main function body.

- End ordinary assignments with semicolons. Short, closely related operations may share a line; do not compress a scientific stage into a dense chain of unrelated statements.

- Wrap long calls and expressions with `...`, aligning continuations with the relevant argument or expression. Align tensor and coordinate rows where that improves readability.

#### 7. End with the scientific result

- Prefer Spinach's plotting conventions: `kfigure`, `plot_1d` or `plot_2d` for spectra, and `kxlabel`, `kylabel`, `klegend`, and `kgrid` where appropriate. Use standard MATLAB plotting when the result is not a standard spectrum.

- Give axes physical labels and units. Make normalisation and comparison conditions explicit; do not silently rescale away a physical difference.

- Use a concise `disp` or `report` for quantities that are better reported numerically. Save a figure, waveform, or data file when it is part of the demonstrated workflow, not as automatic reporting boilerplate.

- Keep assertions and reference comparisons when correctness, convergence, or benchmarking is the subject of the example. Do not graft a generic test harness or progress-reporting framework onto an ordinary simulation.

### Applying the instructions

Choose a nearby example from the same scientific family as the starting point, then apply the rules above. Preserve the existing scientific content when restyling: parameter values, data, references, approximations, ordering, and result interpretation must not change accidentally. Follow the current repository's `AGENTS.md` if its rules change. Legacy deviations are not instructions for new code.

## Documentation and Agent Resources

### Shared Documentation Rules

* **Accurate Description:** Before writing a Wiki narrative, agent skill entry, or knowledge-base entry for a Spinach function, thoroughly analyse its implementation and understand what it does, how it works, and its key algorithms. Only then write a brief but informative explanation of the function’s behaviour and important operational details. Keep Wiki narratives factual and clear; never speculate or introduce information not present in the code.

* **No inconsequential PR or documentation entries:** Do not add mechanically derived or transient metadata to pull-request changes, the shipped knowledge base, the agent skill, or Wiki documentation when it adds no useful explanation and creates avoidable merge conflicts. This includes LOC and source-line counts, indexed-file and aggregate-line counts, generation timestamps, source commit IDs, transient branch/snapshot labels, checkout-specific absolute paths, and duplicate source-header or call-list restatements. Keep stable file paths, signatures, scientific and numerical facts, meaningful references, and substantive descriptions. Commit IDs needed for review or provenance belong in PR discussion or job records, not in shipped documentation. When updating an entry, remove existing low-value metadata without refreshing unrelated content.

### Wiki Instructions

*Spinach* maintains a Wiki for function documentation. When asked (or required) to create or update a function’s Wiki page, apply the shared rules above and the following requirements:

* **Use the Template:** Always use the provided `wiki_template.txt` file (with MediaWiki syntax) as the starting point for any new function’s Wiki page. This template defines the required sections and formatting.

* **Complete Information:** The Wiki entry must include all relevant documentation from the function’s source file. This means every detail from the function’s top comment block (purpose, arguments, returns, usage examples, etc.) must appear on the Wiki page. Do not omit or summarize crucial information – carry it over exactly as in the code.

* **Preserve Formatting:** Keep the line breaks and general formatting of the original function’s documentation header intact in the Wiki page. This ensures consistency between the code and its documentation. For example, if the code’s documentation has separate lines for each parameter, the Wiki should reflect the same line structure.

### Agentic Skill and Knowledge Base

*Spinach* maintains an AI agent skill in `interfaces\agents\spinach-skill` and an AI agent knowledge base in `interfaces\agents\spinach-knowledge`, both to enable third-party AI agents to use Spinach competently. Suggest that the user install both locally. After each Spinach code change, update both appropriately, applying the shared preservation and description rules above. The following additional rules apply to the knowledge base:

* **No line-by-line code regurgitation:** A knowledge entry explains what a function or example does and why; it never restates the source one line at a time. Bullets of the form "Lines 40-41: magnet field; implemented by `sys.magnet=6.9156`" or "Line 18: computes `root_dir` using `root_dir=fileparts(...)`", and sections that list execution stages, control flow, state assignments, or local helper signatures with line numbers, carry no information beyond the source file and are prohibited. Do not add a "Code-derived implementation details" section or any of its subsections to an entry; mechanically generated text of that kind must be deleted before the entry is committed.

* **Local knowledge edits only:** Organise the knowledge base by source-file path: a change to one source file updates its corresponding knowledge entry, if the change warrants documentation, and any other entry whose substantive description is actually affected. Do not ship aggregate category or repository indexes, cross-file registries, generated navigation tables, or other knowledge files that must be refreshed when unrelated source files change. Use the directory layout and file search to discover entries. Do not ship a knowledge entry for the cross-test registry `tests/lib/test_manifest.m`; its membership changes on test additions. Read that source file to inspect registrations. Never recreate the removed section-wide indexes (`etc.md`, `examples.md`, `experiments.md`, `interfaces.md`, `kernel.md`, `tests.md`).

## Confession After Work is Done

After returning the result to the user, you must generate a confession report. Include all aspects of the work that you have skipped, ignored, bypassed, swept under the rug, hallucinated, or did not complete.

