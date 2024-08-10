# [SchrodingerToolkit](https://github.com/walexaindre/SchrodingerToolkit.jl)

This package is designed to solve and develop algorithms for solving certain classes of coupled nonlinear Schrödinger partial differential equations (PDEs). Specifically, we are focused on addressing the following general problem in dimensionless form:

```math
    i \Psi_t = -D \Delta \Psi + \mathcal{V}(\mathbf{x}) \odot \Psi+ \nabla F(|\Psi|^2) \odot \Psi + \mathcal{J}(\Psi)
```

Where: 

* The expected domain of the problem is a bounded Cartesian domain ``\Omega`` from ``\mathbb{R}^d`` with ``d=1,2,3``.

* The boundary conditions are periodic.

* ``\mathbf{x}`` is a vector in the space domain of the problem. Where ``\mathbf{x} := (x_1, x_2, \ldots, x_d)``.

* ``\Psi: \Omega \times (0,T] \longrightarrow \mathbb{C}^N`` with ``\Psi := (\psi_1, \psi_2, \ldots, \psi_N)`` is a complex vector of size ``N``. Here, ``N`` is the number of components of the system with ``N>1``.

* ``|\Psi|^2 := (|\psi_1|^2,|\psi_2|^2,\dots,|\psi_n|^2 )`` represents the vector formed by the squared norm of each component of ``\Psi``.

* ``D:= diag(d_i)`` is a ``\mathbb{R}^{N\times N}`` diagonal matrix with positive dispersive coefficients.

* ``\mathcal{V}:\Omega \longrightarrow \mathbb{R}^{N}`` is an external trapping potential.

* ``\odot`` denotes the Hadamard product.

* ``F:\mathbb{R}^{N} \longrightarrow \mathbb{R}`` is a function which models the strength of intra- and inter-species interactions.

* ``\mathcal{J}:\mathbb{C}^N \longrightarrow \mathbb{C}^N`` models the internal Josephson junction. This Josephson junction takes the form ``\mathcal{J}(\Psi)_n:= \Gamma \displaystyle\sum_{j\neq n}^{N} \psi_j`` where ``\Gamma`` is a real coefficient. 

## What can you expect from this package?

This package is designed to solve the aforementioned problem using a variety of numerical methods. It is modular, allowing users to easily swap out different components of the solution. Additionally, the package is extensible, enabling users to add new numerical methods with ease.

## Important usage considerations

This package is designed as a research tool. As such, it may not be as fast or as stable as other packages. However, it is intended to be flexible and easy to use, allowing users to quickly prototype new algorithms and test new ideas.

## References

For more information about the methods used here you can check the following references:

* [A conservative splitting high-order finite difference method for coupled Gross–Pitaevskii equations in 2D](https://doi.org/10.1140/epjp/s13360-023-04402-6)
* [Finite-difference conservative method for a class of non linear Schrödinger systems](https://doi.org/10.31349/RevMexFisE.19.010205)
* [Finite-difference solutions of a non-linear Schrödinger equation](https://doi.org/10.1016/0021-9991(81)90052-8)
* [Structure preserving Field directional splitting difference methods for nonlinear Schrödinger systems](https://doi.org/10.1016/j.aml.2021.107211)
* [Compact finite difference schemes with spectral-like resolution](https://doi.org/10.1016/0021-9991(92)90324-R)