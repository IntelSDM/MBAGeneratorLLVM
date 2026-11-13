# MBA Equation Generator

This system is designed to simply create unique MBA outputs, not strong ones. We aren't trying to defeat normalization with hardness. <br>
We do this for a VERY VERY important reason, there are 2 ways to defeat MBA resolution. Both of them require equations to be unique, but one requires them to be hard to normalize. <br>
Obviously its all NP, as is any form of maths. It will always be quicker to execute than to normalize any arithmetic equation. <br>
If you are perusing a PHD; trying to find new hardened equations is probably in your interest. <br>
For people who want to do more than mathematics their entire life, the solution to ALL MBA RESOLUTION is simple; opaque predicates. <br>
First of all this is due to truely opaque predicates being impossible to genuinely automatically cancel out, secondly we can generate plausible opaque paths with the truth table. <br>
You need to make opaque predicates actually make sense these days as AI is amazing at resolving repetitive opaque predicates and impossible value brute forcing. This is what this system is for, as you can pollute the path ENTIRELY, while being unresolvable. <br> 
MBA itself is mostly just junkcode these days with the state of normalization optimisation, but it is amazing to be able to use truely possible opaque values in your opaque predicates, to prevent AI based opaque resolution. <br> <br> 
Food for thought for you: most lifting systems operate in passes, over the existing IL, That means in stages.  If you interlayer stages of different systems such as doing an Bitwise opaque predicate mixed with MBA, an optimiser built to handle the MBA to resolve the opaque, interjectively handling them either way round constantly is enough to make lifters drastically slower. Especially with mixed Opaque predicate emission types. 
Think of how the attacker would have to handle some opaques being (MBA, Opaque); (Opaque, MBA, Opaque); (Opaque, MBA); (Opaque type 1, Entirely Different Opaque type 2, MBA), emitting into one another (not exactly ontop but interweaving their requirements on eachother), All in different orders in a single function, your lifter would fail to operate in stages and would have to identify all types manually to optimally attack through optimisation, or the attacker has to brute force passes ontop of eachother with an insane overhead.
<br>That all insists that your opaques aren't runtime dependent; and your battle is between SMT,Z3 optimisation passes from a lifter. The moment runtime opaques enter the image, its different. 

This project is in 3 parts: 
* Const MBA - Bitwise Equations to hide constants.
* Linear MBA - Linear Equations to hide runtime arithmetics.
* Polynomial MBA - Represent ADD,SUB,MUL,SHL constants as a horner encode/decode.

## The Generation logic
* Pick a handful of random boolean expressions: AND/OR/XOR/NOT.
* Evaluate each one at every 0/1 combination of its variables to get a truth table.
* Find coefficients so the weighted sum of the expressions equals the target (a constant or a runtime variable) at every row.
* Solve it with Gaussian elimination. Arithmetic wraps on overflow, so only odd entries can be divided by; any leftover unknowns take random values.
* Plug the coefficients back in. The result is a sum of boolean ops that collapses to the target at every input.

<img width="722" height="941" alt="image" src="https://github.com/user-attachments/assets/83292469-d5e8-4260-b72c-5f6915836c3e" /> <br>
<img width="1094" height="954" alt="image" src="https://github.com/user-attachments/assets/21c4a28a-1301-4576-b525-049512fe568f" /> <br>
<img width="951" height="954" alt="image" src="https://github.com/user-attachments/assets/e48171fe-f956-461b-971a-5fb3d0d99c03" /> <br>
<img width="642" height="910" alt="image" src="https://github.com/user-attachments/assets/27f754e4-dbce-45f2-aa32-15b067b6e92a" /> <br>
<img width="1005" height="651" alt="image" src="https://github.com/user-attachments/assets/cabb1e9e-a194-438f-ac0c-cb98c661b1c4" /> <br>
<img width="716" height="428" alt="image" src="https://github.com/user-attachments/assets/ea33b098-e6df-4565-9491-69beb4c09ee9" /> <br>

## How To Build
* Download [LLVM 19.1](https://github.com/llvm/llvm-project/releases/tag/llvmorg-19.1.0)
* Download [VCPKG](https://github.com/microsoft/vcpkg)
* VCPKG install asmjit
* Open build.bat in a text editor.
* Find the LLVM 19.1 references in the command line, Change them to your file path.
* Find the VCPKG references in the command line, Change them to your file path.
* Run the build script
* Build your application with the ``` -load-pass-plugin ``` Param to the FunctionPass.dll
