# spectral4ode
Spectral Legendre-Galerkin methods for second to fourth order PDE on a rectangular domain with one periodic dimension.

These codes are meant to automate the calculation of neutral stability curves, transition types, and bifucarted solutions. The choice of the type of domain allows us to use a one dimensional spectral method, which makes the execution sufficiently fast to be applied to parametrized equations, i.e. systems where the eigenvalues, eigenfunctions, and bifurcation numbers need to be computed for a large family of different parameter values.

The user needs to implement a function that initializes the problem, initialize_problem.m, which must include the linear parts, parameters, and bilinear interactions. Then a script to run the desired analysis must be written, run_problem.m.

See the examples included for usage instructions on creating initialize_problem.m and run_problem.m files, and the output we can generate.
