# BEPS-SIF-model
The Boreal Ecosystem Productivity Simulator (BEPS) is a two-leaf
enzyme-kinetic model initially developed to estimate the carbon uptake
and the water cycle of boreal ecosystems. It has been substantially upgraded 
since its original release to simulate carbon and water fluxes at an hourly 
time step over various PFTs. Several inter-model comparisons
and site-level validations have shown that BEPS can produce reliable
GPP and ET estimates. 
This version of BEPS model (BEPS-SIF) can simulate the Solar-induced chlorophyll 
fluorescence at the global level. This SIF scheme was developed based on 
radiative  transfer physics to account for the canopy scattering effects 
and can be used at the global scale.

Main modules of BEPS-SIF model:
beps401b_glb_hr: main program for global simulation

inter_prg: This module is a inter-program between main program and modules

evaporation_canopy: This module calculates evaporation and sublimation from 
                    canopy, from overstorey understorey sunlit and shaded leaf

netRadiation: This module calculates net radiation at both canopy level and leaf level

sensible_heat: This module calculates sensible heat from overstorey, understorey and ground

snowpack: This module calculates will calculate sensible heat from overstorey, understorey and ground

surface_temp: This module calculates surface temperature in each step, 
              as well as heat flux for surface to soil layers
              
transpiration: This module calculates transpiration, from overstorey understorey sunlit and shaded leaf

photosyn_gs: This program solves a cubic equation to calculate leaf 
             photosynthesis and chlorophyll fluorescence. Stomatal conductance 
             is computed with the Ball-Berry model. A global sun-induced fluorescence 
             model considers the multiple scattering effects within a canopy.
