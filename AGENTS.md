BACKGROUND INFORMATION FOR ARTIFICIAL INTELLIGENCE AGENTS

==General context==
This is the GitHub repository for Spinach, a spin dynamics simulation software package. Spinach is a scientific computing package that uses the following areas of physics and mathematics: linear algebra, numerical linear algebra, quantum mechanics, Lie algebras, Lie groups, vector and matrix spaces, scientific computing, numerical methods. Spinach covers nuclear magnetic resonance, electron spin resonance, magnetic resonance imaging, quantum optimal control theory, and other areas of science that are related to spin dynamics. Spinach is written exclusively in Matlab.

==Programming style==
Files are structured as single functions with indentation using four spaces. Variables and structure field names use descriptive lowercase names with underscores, for example: zeeman_iso, spin_system, norm_estimate. Every logically distinct code block is preceded by a one-line comment explaining what that code block does.

Documentation is provided as a comment block at the top of each function, explaining the purpose, usage syntax, parameters and outputs. If you are creating a function, you must include this documentation header and follow the same style as you see in other functions residing in "kernel" and "experiments" folders. 

All functions (except those in the example set) call an internal helper named "grumble" right after the main function begins. This helper performs checks on arguments and throws errors if validations fail. If you are creating a function, you must include this validation helper and follow the same style as you see in other functions residing in "kernel" and "experiments" folders.

==Wiki documentation files==
There is a documentation Wiki for this project, it is located at https://spindynamics.org/wiki/ web site. The template for how a Wiki page of each function must look like is in wiki_template.txt file. If you are asked to create a Wiki file for a particular function in Spinach, use that template file. 

When making Wiki documentation file, you must include all description and syntax information that is present in the function itself as well as in its documentaton header. Do not omit any of that information. Keep the line breaks and line folding that you see in the documentation header unchanged when you create a Wiki file.

END OF BACKGROUND INFORMATION
