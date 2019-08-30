// This module will simulate surface temperature in each step, as well as heat flux for surface to soil layers
// edited by XZ Luo, June 1, 2015

/* as it is a interface between ground, air and soil, the core idea is to separate the interface as different layers
by depth of snow, then calculate the temperature gradient and at last calculate the heat flux from ground surface to soil*/

#include "beps.h"
// original beps would use Xg_snow[kkk] at some places after snow melt & frozen, now we uniformly use the value before snow melt & frozen.
void surface_temperature (temp_air,rh_air, depth_snow, depth_water,
                          capacity_heat_soil1, capacity_heat_soil0, Gheat_g, depth_soil1, density_snow,tempL_u,
                          netRad_g, evapo_soil, evapo_water_g, evapo_snow_g,
                          lambda_soil1,percent_snow_g,heat_flux_soil1,
                          temp_ground_last,temp_soil1_last,temp_any0_last, temp_snow_last,
                          temp_soil0_last, temp_snow1_last, temp_snow2_last,
                          temp_ground, temp_any0, temp_snow,
                          temp_soil0,temp_snow1, temp_snow2,
                          heat_flux)

double depth_snow, depth_water;
double density_snow; // density of snow at this step;
double percent_snow_g; // percentage of snow coverage on ground

double lambda_soil1;// thermal conductivity of first layer soil in this step
double depth_soil1;// depth of soil layer 1;
double capacity_heat_soil1,capacity_heat_soil0; // soil heat capacity of first layer soil and 0 layer(?) at step kkk
double heat_flux_soil1; // the heat flux from soil layer 1 to the next layer

double Gheat_g; // aerodynamic conductance of heat at ground, Gheat_g=1/ra_g

double tempL_u;
double netRad_g;// net radiation on grounddouble evapo_soil, evapo_water_g, evapo_snow_g; // evaporation and sublimation from ground
double temp_air, rh_air; // air temperature, relative humidity

double *temp_ground; // ground temperature at this step;
double *temp_any0; // temperature of any layer right above the soil, could be a mixture of snow temperature and soil surface temperature
double *temp_snow; // snow temperature at this step;
double *temp_soil0;// temperature of soil surface right above the soil, the part not covered by snow
double *temp_snow1, *temp_snow2; // temperature of snow layer 2 and 3, used in depth_snow>0.05 m

double temp_ground_last; // ground surface temperature from last step
double temp_any0_last;
double temp_snow_last;
double temp_soil1_last; // temperature of first layer soil in last step
double temp_soil0_last; // temperature of soil surface right above the soil in last step, the part not covered by snow
double temp_snow1_last, temp_snow2_last;

double *heat_flux;//heat flux from ground to soil


{
    double length_step;
    double meteo_pack_output[10];
    double density_air, cp_air; // density of air, specific heat of air
    double cp_ice; // specific heat of ice
    double latent_water, latent_snow;
    double Gg; // radiation available for heating the ground;
    double lambda_snow;// thermal conductivity of snow in this step
    double heat_flux_soil, heat_flux_snow; // heat flux through the soil and snow fraction on ground, separately.
    double heat_flux_snow1, heat_flux_snow2;

    double ra_g; // aerodynamic resistance of heat

    double ttt;// temporary variables

    length_step=kstep;

    cp_ice=2228.261;
    latent_water=(2.501-0.00237*temp_air)*1000000;
    latent_snow=2.83*1000000;
    meteo_pack (temp_air, rh_air, meteo_pack_output);
    density_air = meteo_pack_output[1];
    cp_air = meteo_pack_output[2];
    ra_g=1/Gheat_g;

     /*thermal conductivity of snow*/
    lambda_snow=0.021+4.2*density_snow/10000+2.2*pow(density_snow,3)*pow(10,-9);

    // available energy on ground for
    Gg=netRad_g - evapo_snow_g*latent_snow-(evapo_water_g+evapo_soil)*latent_water;

    //case 1, snow depth<2cm, snow temperature, ground surface temperature, soil surface temperature are the same
    if(depth_snow<=0.02)
	 {
		 ttt=capacity_heat_soil1*0.02/length_step;
		 *temp_ground=(temp_ground_last*(ttt)*ra_g*depth_soil1+Gg*ra_g*depth_soil1+density_air*cp_air*temp_air*depth_soil1+ ra_g*lambda_soil1*temp_soil1_last);
		 *temp_ground=*temp_ground/(density_air*cp_air*depth_soil1+ ra_g*lambda_soil1+ttt*ra_g*depth_soil1);

		 *temp_ground=max(temp_ground_last-25,*temp_ground);
		 *temp_ground=min(temp_ground_last+25,*temp_ground);

		 *temp_any0 = *temp_ground;
		 *temp_snow  =*temp_any0;
		 *temp_soil0 =*temp_any0;
		 *temp_snow1 =*temp_any0;
		 *temp_snow2 =*temp_any0;

		 *heat_flux=2*lambda_soil1*((*temp_any0)-temp_soil1_last)/depth_soil1;
		 *heat_flux=min(100,*heat_flux);
		 *heat_flux=max(-100,*heat_flux);
	 }

    // case 2: depth of snow larger than 2cm, smaller than 5 cm.
    // snow fraction on ground decide the snow temperature based on energy balance
    // soil fraction on ground decide the soil surface temperature based on energy balance
    // snow and soil fraction works in parallel to determine the ground surface temperature
	 else if(depth_snow>0.02  && depth_snow<=0.05)
	 {
		 ttt=capacity_heat_soil1*0.02/length_step; /* for soil fraction part*/

		 *temp_soil0=(temp_soil0_last*ttt*ra_g*depth_soil1+Gg*ra_g*depth_soil1+density_air*cp_air*temp_air*depth_soil1+ 2*ra_g*lambda_soil1*temp_soil1_last)/
			 (density_air*cp_air*depth_soil1+ 2*ra_g*lambda_soil1+ttt*ra_g*depth_soil1);

		 *temp_soil0=max(temp_air-25,(*temp_soil0));
		 *temp_soil0=min(temp_air+25,(*temp_soil0));

		 ttt=cp_ice*density_snow*depth_snow/length_step;  /*for snow part*/
		 *temp_snow=(temp_snow_last*ttt*ra_g*depth_snow+Gg*ra_g*depth_snow+density_air*cp_air*tempL_u*depth_snow+ ra_g*lambda_snow*temp_any0_last)/
			 (density_air*cp_air*depth_snow+ ra_g*lambda_snow+ttt*ra_g*depth_snow);

		 *temp_snow=max(temp_air-25,*temp_snow);
		 *temp_snow=min(temp_air+25,*temp_snow);

		 ttt=(lambda_soil1*temp_soil1_last/depth_soil1+(*temp_snow)*lambda_snow+0.02*capacity_heat_soil1/length_step*temp_any0_last)/
			 (lambda_soil1/depth_soil1+lambda_snow/depth_snow+0.02*capacity_heat_soil1/length_step);
		 *temp_any0=(*temp_soil0)*(1-percent_snow_g)+ttt*percent_snow_g;

		 heat_flux_snow = lambda_snow/(depth_snow+0.5*depth_soil1)*((*temp_snow)-temp_soil1_last);
		 heat_flux_soil = heat_flux_snow*((*temp_any0)-temp_soil1_last)/depth_soil1;

		 *heat_flux=heat_flux_snow*percent_snow_g+heat_flux_soil*(1-percent_snow_g);

		 *heat_flux=min(100,*heat_flux);
		 *heat_flux=max(-100,*heat_flux);

		 /* starting to melt*/
		 if(*temp_snow>zero && temp_snow_last<=zero && depth_snow>zero)
		 {
			 *temp_snow=0;
		 }

		 /* starting to frozen*/
		 if(*temp_snow<zero && temp_snow_last>=zero && depth_water>zero)
		 {
			 *temp_snow=0;
		 }

		 //percent_snow_g =min(1.0,Wg_snow[kkk] / (0.05 * rho_snow[kkk])); // use the fraction before
		 *temp_ground =(*temp_snow)*percent_snow_g+(*temp_soil0)*(1-percent_snow_g);
		 *temp_ground=max(temp_air-25,*temp_ground);
		 *temp_ground=min(temp_air+25,*temp_ground);

		 *temp_snow1=*temp_snow;
		 *temp_snow2=*temp_snow;
	 }

	 // case 3, snow depth is larger than 5 cm
	 // snow coverage on ground is 100%
	 // the first layer of snow is set as 2 cm
	 // the second layer of snow is set as 2 cm, too
	 // the depth of third snow layer is depth_snow-0.04;
	 else if(depth_snow>0.05)
	 {
		 ttt=cp_ice*density_snow*0.02/length_step;

		 *temp_snow=(temp_snow_last*ttt*ra_g*0.04+Gg*ra_g*0.02+density_air*cp_air*temp_air*0.04+ra_g*lambda_snow*temp_snow1_last)/
			 (density_air*cp_air*0.04+ ra_g*lambda_snow+ttt*ra_g*0.04);
		 *temp_snow=max(temp_air-25,*temp_snow);
		 *temp_snow=min(temp_air+25,*temp_snow);

		 heat_flux_snow=lambda_snow*(*temp_snow-temp_snow1_last)/0.04; // why 0.04 here?

		 *heat_flux=heat_flux_snow;
		 *heat_flux=min(100,*heat_flux);
		 *heat_flux=max(-100,*heat_flux);

		 heat_flux_snow1 = lambda_snow*(temp_snow1_last-temp_snow2_last)/(depth_snow-0.02);
		 *temp_snow1=temp_snow1_last+((*heat_flux)- heat_flux_snow1)/(cp_ice*density_snow*0.02 )*length_step ;
		 heat_flux_snow2 = (temp_snow2_last-temp_any0_last)/(0.5*(depth_snow-0.04)/lambda_snow+0.02/lambda_soil1); //??
		 *temp_snow2=temp_snow2_last+(heat_flux_snow1- heat_flux_snow2)/(cp_ice*density_snow*(depth_snow-0.04))*length_step ;

		 *temp_any0 = temp_any0_last+(heat_flux_snow2- heat_flux_soil1) / (capacity_heat_soil0 * 0.02)*length_step;
		 *temp_soil0=*temp_any0;
		 /* starting to melt*/
		 if(*temp_snow>zero && temp_snow_last<=zero && depth_snow>zero)
		 {
			 *temp_snow=0;
		 }

		 /* starting to frozen */
		 if(*temp_snow<zero && temp_snow_last>=zero && depth_water>zero)
		 {
			 *temp_snow=0;
		 }

		 *temp_ground=*temp_snow;
	 }
}

/*==========================*
theoretical top layer: temp_any0
heat_flux_snow2
snow layer 3: temp_snow2
heat_flux_snow1
snow layer 2: temp_snow1
heat_flux_snow
snow layer 1: temp_snow
*/
