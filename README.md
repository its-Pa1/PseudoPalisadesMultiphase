# **Repository Overview – PseudoPalisadesMultiphase**  

Welcome to this GitHub repository, which houses the code associated with the paper **"Multiphase modelling of glioma pseudopalisading under acidosis"**, available at [AIMS Press](https://www.aimspress.com/article/doi/10.3934/mine.2022049), and **Chapter 4** of the thesis titled *Diss_Kumar_Pawan.pdf*, accessible at [KLUEDO](https://kluedo.ub.rptu.de/frontdoor/index/index/docId/6573).  

This repository primarily focuses on the **thesis results**, which provide a **more detailed and comprehensive** presentation of the work compared to the article. The thesis includes all aspects covered in the paper, along with additional insights and extended analyses.  

## **Key Features**  

This repository provides implementations for:  

- **1D and 2D simulations** of tumor microenvironment models  
- **Numerical solutions** of Partial Differential Equations (PDEs)  
- **Generation of plots and videos** for simulation results  

## **Project Structure**  

### **1D_nondim**  
Contains the codes and results corresponding to **Figure 4.10** in the thesis.  

#### **Contents:**  
- `main_nondim.m` – Main script for **1D simulations**  
- `MPM_1D_nondim.avi` and `MPM_1D_nondim.mp4` – Video results  
- `Plots_pattern/` – Directory containing **plot results**  
- `readme.md` – Detailed description of the scripts and methods used  

### **2D_nondim**  
Contains the codes and results for **Figures 4.1 – 4.9** in the thesis.  

#### **Contents:**  
- `main_2D.m` and `main_2D_III.m` – Main scripts for **2D simulations**  
- `compare_chi_three_occ.m`, `compare_K_equal.m`, etc. – Scripts for specific **comparisons and computations**  
- `Plots_diff_chi_double/`, `Plots_diff_chi_zero/`, `Plots_diff_K_equal_I/`, `Plots_diff_K_equal_II/` – Directories containing **plot results**  
- `readme.md` – Detailed description of the scripts and methods used  

### **main_all_chap4.m**  
This script **runs all the main files** required to generate the above-mentioned simulation results (included in **Chapter 4** of the thesis).  

---
Each subdirectory contains a separate README file explaining its contents in detail.

Feel free to explore the directories to dive deeper into specific aspects of the project. If you have any questions or need further clarification, don't hesitate to reach out.  

## **Author**  
**Pawan Kumar**: [@its-Pa1](https://github.com/its-Pa1)  

© 2025, Pawan Kumar. All Rights Reserved.  
