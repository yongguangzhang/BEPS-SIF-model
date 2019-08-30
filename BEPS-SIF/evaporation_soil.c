// this module calculate evaporation from ground surface/top soil, and the evaporation of snow and pond water on surface
// edited by XZ Luo, May 25, 2015

/* input:
air temperature, ground surface temperature, relative humidity of ground (BEPS take it as the air RH)
percentage of snow cover on ground, depth of water, depth of snow
soil water content on first soil layer, porosity of first soil layer
*/

/* output:
evaporation from soil surface
depth of water and snow on ground after evaporation and sublimation
*/
#include "beps.h"
void evaporation_soil(temp_air, temp_g, rh_air, netRad_g, Gheat_g, percent_snow_g,
                       depth_water, depth_snow, mass_water_g, mass_snow_g, density_snow, swc_g, porosity_g,
                       evapo_soil, evapo_water_g, evapo_snow_g)

double temp_air, temp_g, rh_air; // temperature of air, ground and relative humidity of air
double netRad_g; //net radiation on ground
double Gheat_g; // aerodynamic conductance of heat on ground surface

double *percent_snow_g; // percentage of snow on ground
double *depth_water, *depth_snow; // depth of water and snow on ground, after rainfall and snowfall stage 1, before evaporation. output after substracting evapo
double *mass_water_g,*mass_snow_g; // mass of water and snow on ground, output after substracting evapo
double density_snow;//from snowpack stage 1

double swc_g, porosity_g; // soil water content (from last step) and porosity on ground

double *evapo_soil, *evapo_water_g, *evapo_snow_g; // evaporation from soil, evaporation of pond water and snow on surface
//double length_step; // length of step

{
    double meteo_pack_output[10];
    double density_air_g, cp_air_g, vpd_g, slope_vapor_g, psy_air_g; // density of air, specific heat of air, vpd of air, slope of vapor, psychrometer on ground
    double Gwater_g; // conductance of water on soil surface
    double latent_water, latent_snow;
    double density_water;
    double length_step;
	//double change_depth_water,change_depth_snow; // change of water depth in this step

    /////////////////
    meteo_pack (temp_g, rh_air, meteo_pack_output);
    density_air_g = meteo_pack_output[1];
    cp_air_g = meteo_pack_output[2];
    vpd_g = meteo_pack_output[3];
    slope_vapor_g = meteo_pack_output[4];
    psy_air_g = meteo_pack_output[5];

    latent_water=(2.501-0.00237*temp_air)*1000000;
    latent_snow=2.83*1000000;

    density_water=1025.0;
    //////////////////

    Gwater_g=1/(4.0*exp(8.2-4.2*swc_g/porosity_g));
    length_step=kstep;

    // get the percentage of snow
    if (*depth_snow>0.02) *percent_snow_g=1;
    else *percent_snow_g=(*mass_snow_g)/(0.025*density_snow);
    *percent_snow_g=max(*percent_snow_g,0);
    *percent_snow_g=min(*percent_snow_g,1);

    // when there are pond water on ground, there is evaporation from the water
    if (*depth_water>0 && *depth_snow == 0)
	{
        *evapo_water_g=1/(latent_water )*(slope_vapor_g *(netRad_g*0.8-0)+density_air_g *cp_air_g *vpd_g*Gheat_g)/(slope_vapor_g+psy_air_g*(1+(Gheat_g)/0.01));
	}
    else
	{
        *evapo_water_g = 0;
	}
	//*evapo_water_g = max(-0.002/length_step,*evapo_water_g);
	*evapo_water_g = max(0,*evapo_water_g);
	if(*evapo_water_g >0) *evapo_water_g =min(*evapo_water_g,(*depth_water)*density_water/length_step);

	*depth_water=*depth_water-(*evapo_water_g)/density_water*length_step;
    *depth_water=max(0,*depth_water);
    *mass_water_g=*mass_water_g-(*evapo_water_g)*length_step;
//
//    // when there are snow on ground, there is only evaporation from the snow
    if (*depth_snow >0)
	{
        *evapo_snow_g=1/(latent_snow )*(slope_vapor_g *(netRad_g*0.8-0)+density_air_g *cp_air_g *vpd_g*Gheat_g)/(slope_vapor_g+psy_air_g*(1+(Gheat_g)/0.01))*(*percent_snow_g);
	}
    else
	{
        *evapo_snow_g = 0;
	}
	//*evapo_snow_g = max(-0.002/length_step,*evapo_snow_g);
	*evapo_snow_g = max(0,*evapo_snow_g);
	if(*evapo_snow_g >0) *evapo_snow_g =min(*evapo_snow_g,*mass_snow_g/length_step);

    *mass_snow_g=*mass_snow_g-(*evapo_snow_g)*length_step;
    *mass_snow_g=max(*mass_snow_g,0);

    if(*mass_snow_g>0)
    *depth_snow=*depth_snow-(*evapo_snow_g)/density_snow*length_step;
    else
    *depth_snow=0;
//
//
//    // if evaporation from pond and snow have not consumed all energy, then no evaporation from soil surface
//	//change_depth_water=(*evapo_water_g/density_water);
//	//change_depth_snow=(*evapo_snow_g/density_snow);

//
    if ( *depth_water >0 || *depth_snow >0 )
        *evapo_soil=0;
    else
        {
        *evapo_soil=(1-(*percent_snow_g))*1/(latent_water )*(slope_vapor_g*(netRad_g-0)+density_air_g*cp_air_g*vpd_g*Gheat_g)/(slope_vapor_g +psy_air_g*(1 + Gheat_g/Gwater_g));
        *evapo_soil = max(0,*evapo_soil);
        }

}
