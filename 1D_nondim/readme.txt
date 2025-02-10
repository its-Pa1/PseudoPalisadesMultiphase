Main files/scripts:
The file "main_dim.m" is the main matlab script which simulates the 1D form of model 4.4.39 (non-dimensional form). It saves the video of evolution of components and space-time pattern in "Plots_pattern" folder. 



Method: Both diffusion and taxis are discretised in a conservative manner. The first order upwind scheme in conservative form is used for taxis term, while 5 point stencil central difference scheme is applied for diffusion parts. Explicit Euler scheme is used for the time discretisation of both normal and glioma cell equations. IMEX scheme (implicit Euler for diffusion and explicit for source) is used for acidity.


The explicit diffusion  discretisation puts extra stability condition((D*dt/(dx*dx) <0.5)) along with cfl (velocity*dt/dx <1). Incase of any numerical blow-up please try to reduce the time step.

Requirement: Matlab(2020b, it should work for older version also)
