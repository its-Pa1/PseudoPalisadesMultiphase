Main Files/Scripts:

The file `main_i.m` (i = I, III) is the main MATLAB script that calls all the functions and simulates the 2D model (in non-dimensional form, equation 4.4.39 of the thesis) on a square grid and saves the results in the folder `Plots` and `Videos`.

- `main_2D.m`: Simulates three occlusion sites (Initial conditions given by 4.5.1). The corresponding results are included in Figure 4.3 of the thesis.
- `main_2D_III.m`: Simulates two occlusion sites (Initial conditions given by 4.5.2). The corresponding results are included in Figure 4.4 of the thesis.
- `compare_chi_three_occ.m`: Computes the difference between the base case (model 4.4.39) and the one with \(\chi = 0\) (saved in `Plots_diff_chi_zero`) and the difference between the base case model and the one with \(\chi = 2 \times \chi\) (saved in `Plots_diff_chi_double`). The corresponding results are included in Figures 4.5 and 4.6 of the thesis.
- `compare_K_equal.m`: Computes and saves the difference between the base case model 4.4.39 and its counterparts with equal \(K\) (all \(K_{ij}\)) [Experiment 4.4 in Chapter 4 of the thesis]. The results of the three differences (\(K_{ij} = K_{cm}\), \(K_{ij} = K_{mn}\), \(K_{ij} = K_{cn}\)) are saved in `Plots_diff_K_I`, `Plots_diff_K_II`, `Plots_diff_K_III` and `Videos_diff_K_I`, `Videos_diff_K_II`, `Videos_diff_K_III`. The corresponding results are included in Figures 4.7, 4.8, and 4.9 of the thesis.
- `main_all.m`: Calls all the above scripts and saves the 2D results.

Functions:
- `set_const_diff.m`: Assembles the acidity (h) diffusion matrix, implicit part for IMEX scheme.
- `compute_three_occ.m`: Simulates the base model (4.4.39) and saves the involved components.
- `compute_three_occ_chi_zero.m`: Simulates the base model (4.4.39) with \(\chi = 0\) and saves the involved components.
- `compute_three_occ_chi_double.m`: Simulates the base model (4.4.39) with \(2 \times \chi\) and saves the involved components.
- `compute_all_K_equal_I.m`: Simulates the base model (4.4.39) with \(K_{ij} = K_{cm}\) and saves the involved components.
- `compute_all_K_equal_II.m`: Simulates the base model (4.4.39) with \(K_{ij} = K_{mn}\) and saves the involved components.
- `compute_all_K_equal_III.m`: Simulates the base model (4.4.39) with \(K_{ij} = K_{cn}\) and saves the involved components.

Method:
- For time discretization, the explicit Euler method is used for tumor (C) and normal cells (m) equations, and the IMEX method is used for acidity (h). In the IMEX method, the diffusion part is solved implicitly (implicit Euler), while the source terms are solved explicitly (Euler).
- For space discretization, a standard 5-point stencil (central difference) is applied for diffusion. The upwind scheme (first order, conservative form) is used in both x and y directions for the taxis terms (normal and tumor cells).

Note:
The explicit diffusion discretization imposes an additional stability condition \((D \times dt / (dx \times dx) < 0.5)\) along with the CFL condition \((velocity \times dt / dx < 1)\). In case of any numerical blow-up, please try to reduce the time step.

Requirement:
- MATLAB (2020b, it should work for older versions as well)