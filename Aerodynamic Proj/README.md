This project is divided in two main parts:

**Task 1&2**

This folder contains scripts, supporting functions and data files for simulations for task 1 and 2:\
    - the script `Hess_smith_validation.m` uses the file `NACA 0008` for airfoil discretization and compares the Cp distribution obtained from XFOIL (stored in the `Cp` folder) with the results from the Hess-Smith implementation. used outputs: Cl and Cm values, fig. 1 of the report;\
    - the script `cf_comparison.m` compares the skin friction coefficient distribution along the upper surface of two different airfoils (NACA 0008 and NACA 4417) using data stored in `Cf` folder.  used outputs: fig. 2 of the report; 

**Task 1&2 > Bash** folder contains XFOIL bash scripts to obtain comparison data for the previous plots

___________________________________________________________________________________

**Weissinger (Task 3&4)**

!! WARNING: most of these scripts take a significant amount of time to run (several minutes, depending on mesh size/discretization options)

This folder contains scripts, supporting functions and data files for simulations for task 3 and 4:\
       - the script `Task_3_Weissinger_Convergence.m` compares Prandtl’s analytical results with Weissinger’s method results for different system sizes, both plotting and creating a table with the results. used outputs: fig. 3 of the report;\
       - the script  `Task_4_Weissinger_Cessna172_Canadair.m` produces most of the plots seen in the report, comparing the principal aerodynamic characteristics of the wings, tails and whole aircrafts. used outputs: fig. 4-7 of the report.

**Weissinger > Other_Weissinger_Files**

This folder contains scripts which study the impact of isolated geometrical input parameters of an isolated wing while varying the angle of attack. All the files are named according to the studied parameter. These results were not directly included in the report but helped while commenting on the various plots of the Cessna172 and Canadair.

___________________________________________________________________________________


Group Members:

Pasta Emma Maria, 306115\
Rilievi Riccardo, 304220\
Vandoni Giovanni, 307338\
Vescovi Orazio,   308516
