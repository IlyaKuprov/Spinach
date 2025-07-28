INSTRUCTIONS FOR ARTIFICIAL INTELLIGENCE AGENTS


==General context==

This is the GitHub repository for Spinach, a spin dynamics simulation software package. Spinach is a scientific computing package that uses the following areas of physics and mathematics: linear algebra, numerical linear algebra, quantum mechanics, Lie algebras, Lie groups, vector and matrix spaces, scientific computing, numerical methods. Spinach covers nuclear magnetic resonance, electron spin resonance, magnetic resonance imaging, quantum optimal control theory, and other areas of science that are related to spin dynamics. Spinach is written exclusively in Matlab.


==Programming style instructions==

You must create all new functions as standalone files with indentation using four spaces. Use descriptive lowercase names with underscores for all variables and function names, for example: zeeman_isotropic, spin_system, norm_estimate. Above every logically distinct code block, write a one-line comment explaining what that code block does. Do not put comments in the same line after code, comments must always be in the line above the code.

You must provide documentation as a comment block at the top of each function, explaining the purpose, usage syntax, parameters and outputs. If you are creating a function, you must include this documentation header and follow the same style as you see in other functions residing in "kernel" and "experiments" folders. 

You must program all functions (except those in the example set) to call an internal helper named "grumble" right after the main function begins. Put the helper function code at the end of the file. Make the helper function check arguments and throw informative errors if validations fail. If you are creating a new function, you must include this validation helper and follow the same style as you see in helpers of other functions residing in "kernel" and "experiments" folders.

You must leave no spaces around arithmetical ("+","-","*", etc.), logical ("==",">","<=", etc.), or assignment operations anywhere in the code. In general, follow the same style as you see in other functions residing in "kernel" and "experiments" folders. 
 

==Wiki documentation files==

The template for a Wiki page of each function is provided in wiki_template.txt file. If you are asked to create a Wiki file for a function in Spinach, you must use that template file. 

When making a Wiki documentation file for a function, you must include all documentation and syntax information that is present in that function's file. Do not omit any of that information. You must keep the line breaks in the documentation header unchanged when you create a Wiki file. When writing the description of what the function does and how, analyse the function, understnad what it does, and then write a brief informative description. 


==Instruction policies==

You must do everything I ask and produce all the files that I ask you to produce. Keep running until you have produced everything I have told you to produce. Then check my instructions again and make sure you have produced everything I have asked you to produce. If not, then keep running until you you have produced everything I have asked you to produce.



