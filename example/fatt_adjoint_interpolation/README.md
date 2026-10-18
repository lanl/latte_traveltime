# FATT gradients for receivers off the grid

This example compares the FATT gradients of the two settings of `adjoint_receiver` with finite-difference gradients, for receivers between grid points. `nearest` imposes the traveltime residuals as a Dirichlet boundary condition at the grid points nearest to the receivers, and `linear` injects them as true adjoint sources at the receiver locations ([doc/parameters.md](../../doc/parameters.md#dirichlet-condition-or-adjoint-source)).

The model is 1200 m by 600 m on a 10 m grid, with one source at (100, 0) m. One receiver at x = 1000 m moves from z = 200 m to 210 m in 1 m steps, from a grid point through a cell to the next grid point. Two other receivers lie in one cell, at (801.5, 402.0) m and (804.0, 404.5) m, with traveltime residuals of +1 ms and −1 ms. FATT normalizes the model update, so we compare the shape of the gradient near the receivers. The error is $\| \alpha \mathbf{g} - \mathbf{g}^{\mathrm{FD}} \| / \| \mathbf{g}^{\mathrm{FD}} \|$, where $\mathbf{g}$ is the gradient, $\mathbf{g}^{\mathrm{FD}}$ is the finite-difference gradient and $\alpha$ is the least-squares scaling factor.

| Case | `linear` | `nearest` |
|---|---|---|
| Moving receiver on a grid point | 1.8% | 1.8% |
| Moving receiver in the middle of a cell | 1.5% | 64% |
| Receiver pair in one cell | 0.4% | zero gradient |

On grid points, the two settings give the same gradient. Between grid points, `nearest` places the receiver at the nearest grid point, so its error near the receiver grows to 64% in the middle of a cell, whereas `linear` follows the receiver. For the pair, the exact gradient keeps the difference of the two ray paths near the receivers. `nearest` sums the two residuals at their common grid point and gives zero. These differences are local to the receivers, and preconditioning and gradient smoothing reduce them in inversions.

`run.sh` reproduces these results in about one minute and saves `fig_setup.png`, `fig_gradients.png` and `fig_moving_receiver.png`. It needs LATTE binaries in `../../bin`, `numpy` and `matplotlib`.
