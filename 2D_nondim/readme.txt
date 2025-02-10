Main files/scripts:

The file "main_i.m" (i = I,III) is the main matlab script which calls all the functions and simulates the 2-D model (in non-dimensional form, equation 4.4.39 of the thesis) on a square grid and saves the results in the folder Plots and Videos. 

main_2D.m : three occlusion sites (Initial conditions given by 4.5.1), the corresponding results are included in Figure 4.3 of the thesis.

main_2D_III: two occlusion sites (Initial conditions given by 4.5.2), the corresponding results are included in Figure 4.4 of the thesis.

compare_chi_three_occ: difference between the base case (model 4.4.39) and the one with \chi = 0 saved in Plots_diff_chi_zero and the difference between the base case model and the one with \chi = 2*(chi in base case) saved in Plots_diff_chi_double. The corresponding results are included in Figure 4.5 and 4.6 of the thesis.

compare_K_equal: this function computes and saves the difference between the base case model 4.4.39 and the its counterparts with  equal K(all K_ij) [Experiment 4.4 in Chapter 4 of the thesis]. The results of the three differences (K_ij = Kcm, K_ij = Kmn, K_ij = Kcn)  are saved in Plots_diff_K_I, Plots_diff_K_II, Plots_diff_K_III and Videos_diff_K_I, Videos_diff_K_II, Videos_diff_K_III. The corresponding results are included in Figure 4.7, 4.8 and 4.9 of the thesis.

main_all.m: calls all the above scripts and saves the 2D results.




Functions:
set_const_diff.m: assembles the acidity(h) diffusion matrix, implicit part for IMEX scheme.

compute_three_occ: simulates the base model (4.4.39) and saves the involved components
compute_three_occ_chi_zero: simulates the base model(4.4.39) with chi =0 and saves the involved components
compute_three_occ_chi_double: simulates the base model(4.4.39) with 2*chi and saves the involved components
compute_all_K_equal_I: simulates the base model(4.4.39) with K_ij = Kcm and saves the involved components
compute_all_K_equal_II: simulates the base model(4.4.39) with K_ij = Kmn and saves the involved components
compute_all_K_equal_III: simulates the base model(4.4.39) with K_ij = Kcn and saves the involved components

Method: for time discretisation explicit Euler method is used for tumor(C) and normal cells(m) equation and  IMEX method is used for acidity(h). In IMEX method, the diffusion part is solved implicitly (implicit Euler), while the source terms are explicitly (Euler).
 
For space discretisation, standard 5 point stencil( central diff) is applied for diffusion. Upwind scheme (first order, conservative form) is used in both x and y direction for the taxis terms (normal and tumor cells).

The explicit diffusion  discretisation puts extra stability condition((D*dt/(dx*dx) <0.5)) along with cfl (velocity*dt/dx <1). 

Requirement: Matlab(2020b, it should work for older version also)
