Main Files/Scripts:
The file `main_dim.m` is the primary MATLAB script that simulates the 1D form of model 4.4.39 (non-dimensional form). It saves the video of the evolution of components and the space-time pattern in the `Plots_pattern` folder.

Method:
- Both diffusion and taxis are discretized in a conservative manner.
- The first-order upwind scheme in conservative form is used for the taxis term.
- A 5-point stencil central difference scheme is applied for the diffusion parts.
- The explicit Euler scheme is used for the time discretization of both normal and glioma cell equations.
- The IMEX scheme (implicit Euler for diffusion and explicit for source) is used for acidity.

Note:
The explicit diffusion discretization imposes an additional stability condition (D*dt/(dx*dx) < 0.5) along with the CFL condition (velocity*dt/dx < 1). In case of any numerical blow-up, please try to reduce the time step.

Requirement:
- MATLAB (2020b, it should work for older versions as well)