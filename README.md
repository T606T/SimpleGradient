# SimpleGradient Library

SimpleGradient is a lightweight C++ library designed for educational and experimental purposes in numerical optimization. It currently implements a robust one-dimensional gradient descent solver, suitable for understanding core optimization concepts and building a foundation for more advanced methods.

🚧 ## Status: In Development

This library is under active development as part of a personal project to master C++ and numerical optimization. The current version supports single-variable function minimization with derivative-based gradient descent and Armijo backtracking line search.

✨ ## Features

Generic templated solver (SimpleGradient<T>) supporting float and double

Customizable stopping criteria and tolerance levels

Configurable verbosity (INFO, DEBUG, TRACE, ERROR)

Logging support via configurable logger class

Armijo backtracking line search (adaptive step sizing)

Second derivative flatness check (D2Check)

Verbose debug output for educational purposes

🧠 ## Motivation

The goal of this project is to:

Practice modern C++ (RAII, STL, templates, etc.)

Build a reusable and extensible numerical library

Strengthen fundamentals in numerical analysis, optimization, and gradient-based algorithms

Create a visible and credible portfolio piece for R&D and systems programming roles

🛠 Planned Features (Checklist)



📁 ## Structure

SimpleGradient/
├── include/
│   └── Gradient.h      # Header-only implementation of SimpleGradient
├── src/
│   └── test.cpp        # Test runner for predefined functions
├── DebugFile.txt           # Debug log 
├── README.md

🚀 ##Getting Started

$ g++ -Iinclude src/test.cpp -o build/gradient
$ ./build/gradient

📚 ## Example Use

Predefined test cases include quadratic, polynomial, trigonometric, and exponential functions. Each test includes initial guesses and step sizes.

Result Res = Solver.Solve(Functions::Quad, "Quadratic", Derivatives::dQuad, 100, 0.001f, 5, 0.5);

🤝 ## Contributions

Not accepting pull requests yet, but ideas and critiques are welcome.

🧾 ## Author

This project is developed and maintained by a passionate learner aiming to become a C++ developer in R&D-heavy fields such as control systems, embedded development, scientific computing, or HPC.

Feel free to reach out if you're interested in similar topics or have feedback on this project.
