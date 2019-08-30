// this module calculates transpiration, from overstorey understorey sunlit and shaded
// edited by XZ Luo May 20, 2015

/* A transformation of Penmman-Monteith equation is used here.
It could be regarded as a mass transfer process, water vapor inside cells are required by VPD from air and VPD on leaf surface
*/

/* input includes:
temperature of sunlit and shaded leaves from other storey (leaf temperature module);
temperature of air, relative humidity,
conductance of water for sunlit shaded leaves from overstorey and understorey;
leaf area index, sunlit and shaded, overstorey and understorey (from leaf area index module);
*/

/* output:
transpiration from overstorey and understorey
*/
#include "beps.h"
void transpiration (tempL_o_sunlit, tempL_o_shaded, tempL_u_sunlit, tempL_u_shaded, temp_air, rh_air,
                    Gtrans_o_sunlit, Gtrans_o_shaded, Gtrans_u_sunlit, Gtrans_u_shaded,
                    lai_o_sunlit, lai_o_shaded, lai_u_sunlit, lai_u_shaded,
                    trans_o, trans_u)

double tempL_o_sunlit, tempL_o_shaded, tempL_u_sunlit, tempL_u_shaded;
double temp_air, rh_air;
double Gtrans_o_sunlit, Gtrans_o_shaded, Gtrans_u_sunlit, Gtrans_u_shaded;
double lai_o_sunlit, lai_o_shaded, lai_u_sunlit, lai_u_shaded;
double *trans_o, *trans_u;
{
    double LHt_o_sunlit, LHt_o_shaded, LHt_u_sunlit, LHt_u_shaded; // latent heat from leaves W/m2
    double meteo_pack_output[10];
    double density_air, cp_air, vpd_air, slope_vapor_air, psy_air;
    double latent_water;

    meteo_pack (temp_air, rh_air, meteo_pack_output);
    density_air = meteo_pack_output [1];
    cp_air = meteo_pack_output [2];
    vpd_air = meteo_pack_output [3];
    slope_vapor_air = meteo_pack_output [4];
    psy_air = meteo_pack_output [5];

    latent_water=(2.501-0.00237*temp_air)*1000000;

    //////////////////
    LHt_o_sunlit =(vpd_air+slope_vapor_air *(tempL_o_sunlit -temp_air ))*density_air*cp_air*Gtrans_o_sunlit /psy_air;
    LHt_o_shaded =(vpd_air+slope_vapor_air *(tempL_o_shaded -temp_air ))*density_air*cp_air*Gtrans_o_shaded /psy_air;

    LHt_u_sunlit =(vpd_air+slope_vapor_air *(tempL_u_sunlit -temp_air ))*density_air*cp_air*Gtrans_u_sunlit /psy_air;
    LHt_u_shaded =(vpd_air+slope_vapor_air *(tempL_u_shaded -temp_air ))*density_air*cp_air*Gtrans_u_shaded /psy_air;

    //////////////////
    *trans_o =1/(latent_water )*(LHt_o_sunlit *lai_o_sunlit +LHt_o_shaded *lai_o_shaded );
    *trans_u =1/(latent_water )*(LHt_u_sunlit *lai_u_sunlit +LHt_u_shaded *lai_u_shaded );
}
