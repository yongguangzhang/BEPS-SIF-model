// this module will calculate sensible heat from overstorey, understorey and ground
// edited by XZ Luo, May 23, 2015
/* input includes:
temperature of sunlit and shaded leaves from other storey (leaf temperature module);
temperature of air, relative humidity,
temperature of ground (soil heat flux module);
aerodynamic heat conductance of sunlit shaded leaves from overstorey and understorey;
aerodynamic heat conductance of ground;
leaf area index, sunlit and shaded, overstorey and understorey (from leaf area index module);
*/

/* output includes:
sensible heat from overstorey, understorey and ground*/
# include "beps.h"
void sensible_heat (tempL_o_sunlit, tempL_o_shaded, tempL_u_sunlit, tempL_u_shaded, temp_g, temp_air, rh_air,
                    Gheat_o_sunlit, Gheat_o_shaded, Gheat_u_sunlit, Gheat_u_shaded, Gheat_g,
                    lai_o_sunlit, lai_o_shaded, lai_u_sunlit, lai_u_shaded,
                    SH_o, SH_u, SH_g)

double tempL_o_sunlit, tempL_o_shaded, tempL_u_sunlit, tempL_u_shaded; //temperature of overstorey, understorey, sunlit and shaded
double temp_g; // temperature of ground
double temp_air, rh_air;
double Gheat_o_sunlit, Gheat_o_shaded, Gheat_u_sunlit, Gheat_u_shaded; // aerodynamic resistance of heat for each part
double Gheat_g;
double lai_o_sunlit, lai_o_shaded, lai_u_sunlit, lai_u_shaded; // lai
double *SH_o, *SH_u, *SH_g; // output sensible heat
{
    double SH_o_sunlit, SH_o_shaded, SH_u_sunlit, SH_u_shaded; // sensible heat from leaves
    double meteo_pack_output[10];
    double density_air0, cp_air, vpd;

    meteo_pack (temp_air, rh_air, meteo_pack_output);
    density_air0 =meteo_pack_output[1];
    cp_air = meteo_pack_output[2];
    vpd = meteo_pack_output[3];

    //////////////
    SH_o_sunlit=(tempL_o_sunlit-temp_air)*density_air0*cp_air*Gheat_o_sunlit;
    SH_o_shaded=(tempL_o_shaded-temp_air)*density_air0*cp_air*Gheat_o_shaded;

    SH_u_sunlit=(tempL_u_sunlit-temp_air)*density_air0*cp_air*Gheat_u_sunlit;
    SH_u_shaded=(tempL_u_shaded-temp_air)*density_air0*cp_air*Gheat_u_shaded;
    ///////////////

    *SH_o = SH_o_sunlit*lai_o_sunlit+SH_o_shaded*lai_o_shaded;
    *SH_u = SH_u_sunlit*lai_u_sunlit+SH_u_shaded*lai_u_shaded;

    *SH_o=max(-200, *SH_o);
    *SH_u=max(-200, *SH_u);

    *SH_g=(temp_g-temp_air)*density_air0*cp_air*Gheat_g;
}
