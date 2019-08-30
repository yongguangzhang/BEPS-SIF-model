// rainfall_stage1 happens before evaporation of intercepted water from canopy (supply)
// rainfall_stage2 happens after evaporation of intercepted water from canopy (demand)
/*output:
percentage of canopy covered by rainfall, overstorey and understorey (provided to evaporation_canopy);
mass of water available for evaporation on canopy in this step
precipitation on ground
optical output: intercepted mass of rainfall in this step
*/

/*input:
air temperature
percipitation (m/s)
remain of water on leaves from last step (kg/m2) per leaf area
leaf area index of overstorey and understorey, excluding stem.
length of this step (s), if 10min, then it is set as 600
air temperature and humidity
*/

#include "beps.h"
void rainfall_stage1(temp_air, precipitation,mass_water_o_last,mass_water_u_last,
					 lai_o, lai_u, clumping,
					 mass_water_o, mass_water_u,
                     percent_water_o, percent_water_u,precipitation_g)

double temp_air, precipitation; // (C),(m/s)
double mass_water_o_last, mass_water_u_last; // remains of water from last step
double lai_o, lai_u; // leaf area index, no stem
double clumping; //clumping index
//double length_step;
double *mass_water_o, *mass_water_u; // mass of water on leaves, kg/m2, per ground area
double *percent_water_o, *percent_water_u;
double *precipitation_g; // precipitation on ground

{
    double precipitation_o, precipitation_u; // rate of precipitation on overstorey, understorey and ground

    double massMax_water_o, massMax_water_u; // Maximum mass of water could be intercepted per leaf area, in kg/m2
    double massStep_water_o, massStep_water_u; // mass of water intercepted in this step per leaf area, in kg/m2
    double density_water;

    double length_step;

    length_step=kstep;

    if (temp_air < 0|| temp_air == 0) // temperature > 0, otherwise it is snow fall
        precipitation=0;

    density_water=1025.0;
    *precipitation_g=0;

    // overstorey
    precipitation_o=precipitation;
    *mass_water_o=mass_water_o_last+precipitation_o*length_step*density_water*(1-exp(-lai_o*clumping));
    massMax_water_o=0.1*lai_o;

    *mass_water_o=max(0,*mass_water_o);
    *mass_water_o=min(massMax_water_o, *mass_water_o);

    massStep_water_o=(*mass_water_o)-mass_water_o_last;
    massStep_water_o=max(0,massStep_water_o);

    *percent_water_o = (*mass_water_o)/massMax_water_o;
    *percent_water_o = min(1,*percent_water_o);


    // understorey
    precipitation_u=precipitation_o-(massStep_water_o)/density_water/length_step;
    *mass_water_u=mass_water_u_last+precipitation_u*length_step*density_water*(1-exp(-lai_u*clumping));
    massMax_water_u=0.1*lai_u;

    *mass_water_u=max(0,*mass_water_u);
    *mass_water_u=min(massMax_water_u, *mass_water_u);

    massStep_water_u=(*mass_water_u)-mass_water_u_last;
    massStep_water_u=max(0,massStep_water_u);

    *percent_water_u = (*mass_water_u)/massMax_water_u;
    *percent_water_u = min(1,*percent_water_u);


    //ground
    *precipitation_g=precipitation_u-(massStep_water_u)/density_water/length_step;
}

// this module will calculate the water remained on canopy surface after evaporation in this step (used for next step)
// edited by XZ Luo, May 25, 2015

// rainfall_stage1 happens before evaporation of intercepted water from canopy (supply)
// rainfall_stage2 happens after evaporation of intercepted water from canopy (demand)
// Note: rainfall on ground is considered in stage1, and then considered in surface water module (or soil moisture module)

/*input:
mass of water on leaves after precipitation in this step,
evaporation from leaves in this step
*/

/*output:
mass of water on leaves after the evaporation on leaves in this step (this value is transferred to next step)
*/
void rainfall_stage2(evapo_water_o, evapo_water_u, mass_water_o, mass_water_u)

double evapo_water_o, evapo_water_u; // evaporation of intercepted rain in this step, overstorey and understorey, kg/m2/s = mm/s
double *mass_water_o, *mass_water_u; // supply of rain on leaves, already added precipitation in this step

{
    double length_step; // length of step
    length_step=kstep;

    *mass_water_o=*mass_water_o-evapo_water_o*length_step;
    *mass_water_o=max(0,*mass_water_o);

    *mass_water_u=*mass_water_u-evapo_water_u*length_step;
    *mass_water_u=max(0,*mass_water_u);
}

