# Ocean-Nitrogen-Isotopes-Box-Models

Prognostic multi-box (5 and 9) ocean model to constrain the systematics of N isotopes at global scale 

The model is run on MATLAB. 

Reference:
Fripiat, F., D.M. Sigman, A. Martínez-García, D. Marconi, X.E. Ai, A. Auderset, S.E. Fawcett, S. Moretti, A.S. Studer and G.H. Haug (2023). The impact of incomplete nutrient consumption in the Southern Ocean on global mean ocean nitrate δ15N. Accepted in Global Biogeochemical Cycles. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FILE DESCRIPTION for the five-box model:

- fivebox_ocean_isotopes_main_template.m : master file for the five-box model which is run with one given set of parameters.
- fivebox_ocean_isotopes_ode.m : ODE solver for the "fivebox_ocean_isotopes_main_template.m"
- fivebox_ocean_isotopes_main_template_sensitivity.m : master file for the sensitivity experiments for the five-box model 
- fivebox_ocean_isotopes_sensitivity_ode.m : ODE solver for the "fivebox_ocean_isotopes_main_template_sensitivity.m"
- fivebox_sensitivity_generation.m : script to generate the "sensitivity" matrix being used in "fivebox_ocean_isotopes_main_template_sensitivity.m"


HOW TO for the five-box model which is run with one given set of parameters:
- to set the model parameters in "fivebox_ocean_isotopes_ode.m"
- to run "fivebox_ocean_isotopes_main_template.m"

HOW TO for the five-box model with the sensitivity experiments:
- to run "fivebox_sensitivity_generation.m" to generate the "sensitivity" matrix by randomly varying parameters over a range well beyond literature estimates 
(see table 1 in the reference manuscript)
- to run "fivebox_ocean_isotopes_main_template_sensitivity.m"

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FILE DESCRIPTION for the nine-box model:

- ninebox_ocean_isotopes_main_template.m : master file for the nine-box model which is run with one given set of parameters.
- ninebox_ocean_isotopes_ode.m : ODE solver for the "ninebox_ocean_isotopes_main_template.m"
- ninebox_ocean_isotopes_main_template_sensitivity.m : master file for the sensitivity experiments for the nine-box model 
- ninebox_ocean_isotopes_sensitivity_ode.m : ODE solver for the "ninebox_ocean_isotopes_main_template_sensitivity.m"
- ninebox_sensitivity_generation.m : script to generate the "sensitivity" matrix being used in "ninebox_ocean_isotopes_main_template_sensitivity.m"


HOW TO for the nine-box model which is run with one given set of parameters:
- to set the model parameters in "ninebox_ocean_isotopes_ode.m"
- to run "ninebox_ocean_isotopes_main_template.m"

HOW TO for the nine-box model with the sensitivity experiments:
- to run "ninebox_sensitivity_generation.m" to generate the "sensitivity" matrix by randomly varying parameters over a range well beyond literature estimates 
(see table 1 in the reference manuscript)
- to run "ninebox_ocean_isotopes_main_template_sensitivity.m"
