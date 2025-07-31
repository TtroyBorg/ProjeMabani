##  Problem Setting

We aim to minimize a real-valued function:

\[
f: \mathbb{R}^n \rightarrow \mathbb{R}, \quad \min_{\mathbf{x} \in \mathbb{R}^n} f(\mathbf{x})
\]

using an iterative method based on **BFGS**, a quasi-Newton technique that builds an approximation to the inverse Hessian.

---

##  BFGS Algorithm Overview

At each iteration:

1. Compute the descent direction:

\[
\mathbf{p}_k = -B_k^{-1} \nabla f(\mathbf{x}_k)
\]

2. Perform a **line search** to find a step size \( \alpha_k \) satisfying the **Goldstein conditions**.

3. Update position:

\[
\mathbf{x}_{k+1} = \mathbf{x}_k + \alpha_k \mathbf{p}_k
\]

4. Update the inverse Hessian approximation \( B_k \) via:

\[
\mathbf{s}_k = \mathbf{x}_{k+1} - \mathbf{x}_k, \quad
\mathbf{y}_k = \nabla f(\mathbf{x}_{k+1}) - \nabla f(\mathbf{x}_k)
\]

\[
B_{k+1} = B_k + \frac{\mathbf{y}_k \mathbf{y}_k^\top}{\mathbf{y}_k^\top \mathbf{s}_k} - \frac{B_k \mathbf{s}_k \mathbf{s}_k^\top B_k^\top}{\mathbf{s}_k^\top B_k \mathbf{s}_k}
\]

---

##  Goldstein Line Search

We select a step size \( \alpha \) such that:

\[
f(\mathbf{x}_k + \alpha \mathbf{p}_k) \le f(\mathbf{x}_k) + \rho \alpha \nabla f(\mathbf{x}_k)^\top \mathbf{p}_k \tag{(1)}
\]
\[
f(\mathbf{x}_k + \alpha \mathbf{p}_k) \ge f(\mathbf{x}_k) + (1 - \rho) \alpha \nabla f(\mathbf{x}_k)^\top \mathbf{p}_k \tag{(2)}
\]

for some \( 0 < \rho < 0.5 \). These ensure the step is **not too short** and **not too long**.

---

##  References

Nocedal & Wright, Numerical Optimization
