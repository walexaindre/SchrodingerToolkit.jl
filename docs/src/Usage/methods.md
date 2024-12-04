# How to choose a method

Since the process of method generation and selection is not fully automated at this time, you need to choose a method based on the specifics of your problem.

Here is a description of available methods and characteristics of each method.

| Method Name | PDE Type | Method Type            | Linear Solver |
|-------------|----------|------------------------|---------------|
| M1          | CNLS     | Crank Nicolson        | GMRES         |
| M2          | CGPE     | Directional Decomposition | GMRES         |
| M3          | CGPE     | Directional Decomposition | LU            |
| M4          | CNLS     | Crank Nicolson        | LU            |