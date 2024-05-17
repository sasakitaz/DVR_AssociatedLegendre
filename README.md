# What is this?
Discrete variable representation (DVR) calculation of spherical harmonics basis

This program solve time independent Schrodinger equation
```math
\begin{eqnarray}
  H = B {\bf{j}}^2 + V(\cos \theta)
\end{eqnarray}
```
by $\ket{j, m}$ basis.

# Output
1. dimension: Hamiltonian matrix size
2. l_potential: Legendre quantum number in potential function of analytical solution
3. grid: numpy library, np.polynomial.legendre.leggauss() and calculated grid value (latter includes former's values)
4. analytical solution by spherical harmonics basis: eigenvalue of Hamiltonian which has analytical matirix element, $V(\cos\theta) = P_{l_{\mathrm{potential}}}(\cos\theta)$
5. numerical solution by associated Legendre-DVR: eigenvalue of Hamiltonian which has DVR matirix element, $V(\cos\theta) = \mathrm{potential function}(\cos\theta)$
