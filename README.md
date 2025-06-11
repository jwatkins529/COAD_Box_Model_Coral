# COAD_Box_Model_Coral

## Description

This repository contains MATLAB scripts to compute the &delta;<sup>13</sup>C, &delta;<sup>18</sup>O, Δ<sub>47</sub>, and Δ<sub>48</sub> of coral aragonite. The full description of the model and parameters can be found in: Watkins, J. Jia, Q., Zhang, S., Devriendt, L., and Chen, S., 2025, A coral biomineralization model for dual clumped isotopes, submitted to Geochemistry, Geophysics, Geosystems. 

![Figure_1](https://github.com/user-attachments/assets/0c264f2d-3dd5-48f6-9728-b7bcc3773c14)


## Requirements
In order to run the scripts you will need Matlab. The codes were produced using Matlab2022a but earlier and later versions will probably work just as well. 

## Scripts

### 1. COAD_Box_Model_Coral.m

The COAD box model allows one to calculate kinetic isotopes effects in the full CaCO<sub>3</sub>-DIC-H<sub>2</sub>O system. The input parameters can be changed in the top section of COAD_Box_Model_Coral.m. The script Run_COAD_Box_Model_Coral.m executes COAD_Box_Model_Coral.m for a range of F<sub>Alk</sub> values and was used to produce the modeling results of the paper. The parameters are currently configured to output the best-fitting curve for deep-sea corals. 
<ul>
  <li> CaCO3_DIC_Coral.m - this is the ion-by-ion function script
  <li>COAD_Box_Model_Coral.m - this is the main script that includes all of the inputs and uses the ODE solver.
  <li>Run_COAD_Box_Model_Coral.m - this is used to run COAD_Box_Model_Coral.m and post-process the results
  <li>Davies_cold.txt - data for deep-sea corals from Davies et al. (2022)
  <li>Davies_warm.txt - data for tropical corals from Davies et al. (2022)

</ul>  
To run, put all files in same directory and execute Run_COAD_Box_Model_Coral.m. This will run the code to steady state for a range of F<sub>Alk</sub> values that are specified on line 4 of Run_COAD_Box_Model_Coral.m.  The output is two figures: (1) a 12-panel figure showing how composition of the calcifying fluid varies with F<sub>Alk</sub>  and (2) a 2-panel figure showing the &delta;<sup>13</sup>C-&delta<sup>18</sup>O and Δ<sub>47</sub>-Δ<sub>48</sub> covariations.

![COAD_output_1](https://github.com/user-attachments/assets/3ef12242-5ab6-4d1c-86f3-d329f26faec2)

![COAD_output_2](https://github.com/user-attachments/assets/8e80fe51-bf72-4ce8-885f-bed963857dbe)

For the paper, we color-coded the curves by parameter (e.g., pH) by saving the outputs to a .txt file and using a different software (GMT).  

![Figure_9](https://github.com/user-attachments/assets/1c764b28-5aba-429b-837d-986fb9717e4c)
