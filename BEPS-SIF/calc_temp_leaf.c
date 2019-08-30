/*************************************************************************
* Subroutine to calculate the sunlit and shaded leaf temperatures for 
* overstory and understory leave. 
*
* Written and refactored by Liming He (liming.he@gmail.com)
* Original contributor: Weimin Ju
* 
* Last update: 09/15/2015
* This subroutine was created on May 15, 2015
**************************************************************************/

#include "beps.h"

void Leaf_Temperatures(double Tair, double slope, double psychrometer, double VPD_air, double Cp_ca, \
	double Gw_o_sunlit, double Gw_o_shaded, double Gw_u_sunlit, double Gw_u_shaded, \
	double Gww_o_sunlit, double Gww_o_shaded, double Gww_u_sunlit, double Gww_u_shaded, \
	double Gh_o_sunlit, double Gh_o_shaded, double Gh_u_sunlit, double Gh_u_shaded, \
	double Xcs_o, double Xcl_o, double Xcs_u, double Xcl_u, \
	double radiation_o_sun, double radiation_o_shaded, double radiation_u_sun, double radiation_u_shaded, \
	double *Tc_o_sunlit, double *Tc_o_shaded, double *Tc_u_sunlit, double *Tc_u_shaded)
{
	*Tc_o_sunlit = Leaf_Temperature(Tair, slope, psychrometer, VPD_air, Cp_ca, \
		Gw_o_sunlit, Gww_o_sunlit, Gh_o_sunlit, Xcs_o, Xcl_o, radiation_o_sun);

	*Tc_o_shaded = Leaf_Temperature(Tair, slope, psychrometer, VPD_air, Cp_ca, \
		Gw_o_shaded, Gww_o_shaded, Gh_o_shaded, Xcs_o, Xcl_o, radiation_o_shaded);

	*Tc_u_sunlit = Leaf_Temperature(Tair, slope, psychrometer, VPD_air, Cp_ca, \
		Gw_u_sunlit, Gww_u_sunlit, Gh_u_sunlit, Xcs_u, Xcl_u, radiation_u_sun);

	*Tc_u_shaded = Leaf_Temperature(Tair, slope, psychrometer, VPD_air, Cp_ca, \
		Gw_u_shaded, Gww_u_shaded, Gh_u_shaded, Xcs_u, Xcl_u, radiation_u_shaded);
}



/* subroutine to calculate leaf temperature */
double Leaf_Temperature(double Tair, double slope, double psychrometer, double VPD_air, double Cp_ca, \
	double Gw, double Gww, double Gh, double Xcs, double Xcl, double radiation)
{
	double p_star;
	double Tc;

	double R = 1 / Gw + 1 / (Gww* (Xcs + Xcl)); 
	p_star = (Gw + Gww* (Xcs + Xcl) ) / psychrometer; 

	Tc = Tair + (radiation - VPD_air * rho_a * Cp_ca * p_star) / (rho_a * Cp_ca * (Gh + slope * p_star) ); //temperature of sunlit leaves, overstorey

	Tc = max(Tair - 3.0, Tc); 
	Tc = min(Tair + 5.0, Tc);

	return Tc;
}
