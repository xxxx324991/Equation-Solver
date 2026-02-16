Equation-Solver (C++)

A robust numerical engine designed to solve quadratic, cubic, and quartic equations. The core of this project is a comprehensive, custom-built Complex Number class that handles the intricate logic of algebraic field extensions.

Key Achievement: Custom Complex Class

While many use standard libraries, this project implements a full-featured Complex Number system from the ground up to ensure maximum control over precision and logic flow.


Technical Highlights:

Comprehensive Operator Overloading: Implemented +, -, *, /, +=, -=, and << for seamless arithmetic.

Precision Handling: Optimized for complex coefficients where standard numerical solvers often struggle.

Robust Logic: Beyond classroom exercises, this class manages real-imaginary interactions with high mathematical rigor.


🚀 Solvability Scope

Quadratic Equations: $ax^2 + bx + c = 0$

Cubic Equations: $ax^3 + bx^2 + cx + d = 0$

Quartic Equations: $ax^4 + bx^3 + cx^2 + dx + e = 0$ (Implemented using radical-based numerical methods).


💻 Code Snippet: Operator Overloading

Here is a glimpse of the logical depth in the Complex class:

Complex operator*(const Complex& other) const
