# Anisotropic Diffusion in 2D

## Problem description

We consider the steady anisotropic diffusion equation on the unit square $(0,1)^2$:

```math
    -u_{xx} - \varepsilon\,u_{yy} = f(x,y),
```

with Dirichlet boundary conditions:

```math
\begin{cases}
u(0,y) = u(1,y) = 4, \\[4pt]
u(x,0) = \sin(\pi x), \\[4pt]
u(x,1) = \alpha\,\sin(\pi x),
\end{cases}
\qquad \text{with } \alpha = 0.8.
```

The source term is defined as

```math
f(x,y) = b_0 \; \sin(\pi x)\sin(\pi y), \qquad b_0 = 1.2.
```

The anisotropy parameter, $\varepsilon$, controls the relative strength of diffusion in the $y$-direction compared to $x$. The values of $\varepsilon$ range from $1$ (isotropic case) down to $10^{-6}$ (highly anisotropic).

---

## Analytical representation

No closed-form elementary solution exists for this mixed Dirichlet problem.  
However, the analytical solution can be expressed as a **Fourier–hyperbolic series**:

```math
\begin{aligned}
u(x,y) = 4 &+ \frac{b_0}{\pi^2(1+\varepsilon)}\,\sin(\pi x)\sin(\pi y) \\
&+ \sum_{n=1}^{\infty}
\sin(n\pi x)\,
\frac{ a_n^{(0)}\,\sinh\!\big(\tfrac{n\pi}{\sqrt{\varepsilon}}(1-y)\big)
     + a_n^{(1)}\,\sinh\!\big(\tfrac{n\pi}{\sqrt{\varepsilon}}y\big)}
     {\sinh\!\big(\tfrac{n\pi}{\sqrt{\varepsilon}}\big)},
\end{aligned}
```

where the boundary coefficients are

```math
a_n^{(0)} = 2\!\int_0^1 \!\!\big(\sin(\pi x)-4\big)\sin(n\pi x)\,dx, \qquad
a_n^{(1)} = 2\!\int_0^1 \!\!\big(\alpha\sin(\pi x)-4\big)\sin(n\pi x)\,dx.
```

This representation is used as a **validation reference** through numerical convergence rather than symbolic evaluation.

---

## Discretisation

The PDE is discretized using a **five-point finite difference stencil** on a uniform grid of size  $N\times N$ with spacing $h = 1/(N+1)$:

```math
\begin{aligned}
&\frac{-u_{i-1,j} + 2u_{i,j} - u_{i+1,j}}{h^2}
+
\varepsilon\,
\frac{-u_{i,j-1} + 2u_{i,j} - u_{i,j+1}}{h^2}
= f_{i,j}.
\end{aligned}
```

The resulting linear system $A \cdot u = b$ is symmetric and positive definite (SPD). The coefficient matrix $A$ exhibits strong ill-conditioning for small $\varepsilon$.

---

## Solvers implemented

<div align="center">
  <table>
    <thead>
      <tr>
        <th align="center">Method</th>
        <th align="center">Description</th>
      </tr>
    </thead>
    <tbody>
      <tr>
        <td align="center"><b>CG</b></td>
        <td align="center">Standard Conjugate Gradient for SPD systems</td>
      </tr>
      <tr>
        <td align="center"><b>PCG (IC(0))</b></td>
        <td align="center">Preconditioned CG with Incomplete Cholesky (IC(0)) factorization</td>
      </tr>
    </tbody>
  </table>
</div>

## Experiment setup

Each solver was tested for interior grid sizes

```math
N \in \{64,\,128,\,256,\,512\},
```

and anisotropy levels

```math
\varepsilon \in \{1, 7\times10^{-1}, 5\times10^{-1}, 3\times10^{-1}, 10^{-1}, 7\times10^{-2}, 5\times10^{-2}, 3\times10^{-2}, 10^{-2}, \dots, 10^{-6}\}.
```

For each configuration, we log:

- iteration count to convergence ($10^{-6}$ tolerance),
- runtime in milliseconds,
- estimated condition number $\hat{\kappa} = \lambda_{\max}/\lambda_{\min}$ via a Lanczos procedure,
- residual norm at termination.

## Results

### Iterations vs $\varepsilon$

<p align="center">
  <img src="./docs/iters_eps_64.png" width="45%">
  <img src="./docs/iters_eps_128.png" width="45%"><br>
  <img src="./docs/iters_eps_256.png" width="45%">
  <img src="./docs/iters_eps_512.png" width="45%">
</p>

### Runtime vs $\varepsilon$

<p align="center">
  <img src="./docs/time_eps_64.png" width="45%">
  <img src="./docs/time_eps_128.png" width="45%"><br>
  <img src="./docs/time_eps_256.png" width="45%">
  <img src="./docs/time_eps_512.png" width="45%">
</p>

### Condition number estimate

<p align="center">
  <img src="./docs/kappa_eps_CG.png" width="60%">
</p>

### Effect of $\varepsilon$ to the solution $u(x, y)$

<p align="center">
  <img src="./docs/heatmap_1e-3.png" width="30%">
  <img src="./docs/heatmap_2e-2.png" width="30%">
  <img src="./docs/heatmap_4e-2.png" width="30%">
  <img src="./docs/heatmap_2e-1.png" width="30%">
  <img src="./docs/heatmap_1e-1.png" width="30%">
  <img src="./docs/heatmap_1.png" width="30%">
</p>

---

## Key observations

- As $\varepsilon \downarrow 0$, diffusion along $y$ weakens, producing **anisotropic stretching** of the solution contours and **increasing ill-conditioning** of the linear system.
- The **condition number** $\hat{\kappa}$ grows roughly as $O(1/\varepsilon)$, causing a sharp rise in CG iteration counts.
- The **IC(0) preconditioner** effectively limits this growth, cutting iteration counts by factors of 2–5 for moderate anisotropy levels.
- For extreme anisotropy ($\varepsilon < 10^{-4}$), the system becomes nearly singular in the $y$-direction, and performance saturates.

---

**Author:** Gabriela Raluca Stancu  
University of Western Macedonia  
Kozani, Greece (2025)
