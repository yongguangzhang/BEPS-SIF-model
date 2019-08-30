/*************************************************************************
* Subroutine to calculate the Vcmax and Jmax for sunlit and shaded big-leaf
*
* Written and refactored by Liming He (liming.he@gmail.com)
* Original contributor: Gang Mo
* Reference: 
*	(1) Chen, J. M., G. Mo, J. Pisek, F. Deng, M. Ishozawa, D. Chan, 2012. 
*	Effects of foliage clumping on global terrestrial gross primary 
*	productivity. Global Biogeochemical Cycles, VOL. 26, GB1019, 18, 
*	doi:10.1029/2010GB003996
*	(2) Medlyn, B.E. et al., 1999. Effects of elevated [CO2] on photosynthesis 
*   in European forest species: a meta-analysis of model parameters. 
*   Plant, Cell & Environment, 22(12): 1475-1495.
*
* Last update: May 5, 2015
* This subroutine was created on Mar. 18, 2014
**************************************************************************/

#include "beps.h"

void Vcmax_Jmax(double lai_o, double clumping, double Vcmax0, \
	double slope_Vcmax_N, double leaf_N, double CosZs, \
	double *Vcmax_sunlit, double *Vcmax_shaded, double *Jmax_sunlit, double *Jmax_shaded)
{
	// Note: Vcmax0 is for the leaf vcmax at top of the canopy. LHE

	double Kn = 0.3;  // Kn = 0.713/2.4;
	double G_theta = 0.5;
	double K,expr1,expr2,expr3;

	if(lai_o < 0.001) lai_o = 0.001; // to avoid error in the code, a minimum is setup here. LHE.

	if (CosZs>0) 
	{
		K = G_theta*clumping/CosZs; // always > 0 LHE

		expr1 = 1 - exp(-K*lai_o); // always > 0. LHE
		expr2 = 1 - exp(-(Kn + K)*lai_o);
		expr3 = 1 - exp(-Kn*lai_o);

		if(expr1>0) 
		  *Vcmax_sunlit = Vcmax0 * slope_Vcmax_N * leaf_N * K*expr2 / (Kn + K) / expr1; // LHE. May 15, 2014@PGB205b
		
		else *Vcmax_sunlit = Vcmax0;

		if (K > 0 && lai_o > expr1/K) 
			//Vcmax_shaded = Vcmax0*parameter1[47]*parameter1[46]*(expr3/Kn-expr2*clumping/(Kn+K))/(lai_o-expr1*2*CosZs); /* modified by Ting, Mar,28,2013 */
			
			*Vcmax_shaded = Vcmax0 * slope_Vcmax_N * leaf_N * (expr3 / Kn - expr2 / (Kn + K)) / (lai_o - expr1/K); // original code.
		//	*Vcmax_shaded = Vcmax0 * slope_Vcmax_N * leaf_N * (expr3 / Kn - expr2 * clumping / (Kn + K)) / (lai_o - 2 * CosZs * expr1); // LHE May 15, 2014@PGB205b.
		        
		else
			*Vcmax_shaded = Vcmax0;
			
	}
	else
	{
		*Vcmax_sunlit = Vcmax0;
		*Vcmax_shaded = Vcmax0;
	}

	*Jmax_sunlit = *Vcmax_sunlit * 2.39 - 14.2;
	*Jmax_shaded = *Vcmax_shaded * 2.39 - 14.2;
}
