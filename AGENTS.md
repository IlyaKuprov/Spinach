# Instructions for Artificial Intelligence Agents

## General Context

*Spinach* is an open-source spin dynamics simulation library implemented exclusively in **MATLAB**. It spans numerous areas of physics and mathematics, including linear algebra (and numerical linear algebra), quantum mechanics, Lie algebras and Lie groups, as well as general scientific computing and numerical methods. *Spinach* supports applications such as nuclear magnetic resonance (**NMR**), electron spin (paramagnetic) resonance (**ESR/EPR**), magnetic resonance imaging (**MRI**), quantum optimal control theory, and other spin dynamics-related domains. This repository contains the *Spinach* codebase, and all contributions or AI-generated code **must** adhere to the established conventions of this codebase.

## Programming Style Guidelines

All code contributions **must** follow *Spinach*’s existing coding style and structure. When writing code, adhere to the following rules without exception:

* **Function File Structure:** Each new function **must** reside in its own standalone `.m` file. Use **four spaces** for indentation (no tabs).

* **Naming Conventions:** Use descriptive, all-lowercase names with underscores for **all** variables and function names. For example, follow naming patterns seen in the codebase such as `zeeman_isotropic`, `spin_system`, or `norm_estimate`. **Never** use capital letters or ambiguous names.

* **Code Comments:** Above every conceptually distinct block of code, write a one-line comment explaining its purpose. **Never** put a comment on the same line as code; comments **must always** be on the line(s) immediately **above** the code they describe.

* **Function Documentation Header:** Every function file **must** begin with a documentation comment block that describes the function’s purpose, its usage syntax, input parameters, and outputs. **Always** format this documentation header exactly as seen in existing functions (refer to the `kernel` and `experiments` directories for examples). Do not omit any expected sections in the header.

* **Input Validation with `grumble`:** All non-example functions **must** perform input argument validation at the start of the function using the `grumble` helper. Immediately after the function definition, call an internal helper function named `grumble` to check the validity of arguments. Define the full `grumble` helper at the end of the same file.

* **Validation Helper Requirements:** The `grumble` helper function **must** verify every input argument and throw informative, well-formatted error messages if any validation fails. Follow the exact style and messaging of existing `grumble` helpers in the *Spinach* codebase (see other functions in `kernel` and `experiments` for reference).

* **Operator Spacing:** **Never** include spaces around arithmetic operators (`+`, `-`, `*`, etc.), logical operators (`==`, `>`, `<=`, etc.), or the assignment operator (`=`). Write expressions like `a=b+c*d` **without** spaces. This convention is consistent across the entire codebase.

* **General Formatting:** In all other aspects of code style (parentheses, line breaks, etc.), **mimic the existing code**. Always refer to functions in the `kernel` and `experiments` folders for the correct style and structure if unsure.

* **Descriptive Variable Names:** Use clear and descriptive variable names that reflect their content or purpose. **Do not** use vague names. The only exceptions are simple loop indices (e.g., `n`, `k` for loop counters). **Do not** use `i` and `l` as variables. 

* **Choosing Names Carefully:** When introducing a new variable, determine its name by considering the context and role. Read the function’s documentation and understand what the function does and what the variable represents. **Then** choose a concise name that conveys that meaning.

* **Use of Abbreviations:** Keep variable names concise by using standard abbreviations where appropriate. For example, a variable holding a **property index** may be named `prop_idx`. Ensure any abbreviation used is commonly understood or documented in the codebase.

## Wiki Documentation Instructions

*Spinach* maintains a Wiki for function documentation. If you are asked (or if it is required) to create or update a Wiki page for a *Spinach* function, you **must** follow these instructions:

* **Use the Template:** Always use the provided `wiki_template.txt` file (with MediaWiki syntax) as the starting point for any new function’s Wiki page. This template defines the required sections and formatting.

* **Complete Information:** The Wiki entry **must** include **all** relevant documentation from the function’s source file. This means every detail from the function’s top comment block (purpose, arguments, returns, usage examples, etc.) **must** appear on the Wiki page. **Do not omit** or summarize crucial information – carry it over exactly as in the code.

* **Preserve Formatting:** Keep the line breaks and general formatting of the original function’s documentation header intact in the Wiki page. This ensures consistency between the code and its documentation. For example, if the code’s documentation has separate lines for each parameter, the Wiki should reflect the same line structure.

* **Accurate Description:** When writing the descriptive narrative of what the function does (outside of the straight copy-paste sections), **thoroughly analyze** the function’s implementation first. Make sure you understand how it works and what its key algorithms are. **Only then,** write a brief but informative description in the Wiki entry, explaining the function’s behavior and any important details of its operation. Always maintain a factual and clear tone – **never** speculate or introduce information not present in the code.

## Execution Policies

* **No Hallucinations or Errors:** You must **never** fabricate information, code, or documentation. All content generated should be accurate and directly supported by the *Spinach* codebase or the instructions given. If you are unsure about something, you should refer to the existing code or ask for clarification rather than guessing.

* **Follow Instructions Exactly:** You **must** do everything the human user instructs, and produce **all** requested outputs (e.g. multiple functions or files) exactly as specified. **Never** ignore any part of the request. Each task in the instructions should be completed fully.

* **Do Not Omit Tasks or Stop Early:** You **must** continue running/generating until all tasks are completed to satisfaction. You **must not** terminate or cut off the output prematurely. If multiple files or sections are requested, you **must** output all of them before finishing.

* **Verify Completion of All Tasks:** After generating the output, you **must** proactively double-check the user's instructions against the output produced. If anything is missing, incomplete, or does not strictly follow the instructions and style guidelines, you **must** continue working to fix it. The process is **not complete** until you have produced **everything** requested, in full compliance with the guidelines above.

