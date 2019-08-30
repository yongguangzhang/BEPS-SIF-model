/*********************************************************************
----------------------------
Program: SOIL_WATER_FACTOR.C
Version 2.0
----------------------------
Description:  Compute soil water stress factor. please refer to file soil_water_factor.cpp for
the orginal version. LHE.

Last revision by LHE. May 22, 2015. 
Rewritten by: Liming He. Jan 29, 2013.
Modified by: Mustapha El Maayar - March 2008
Written by: Weiming Ju
*********************************************************************/
# include "stdio.h"
# include "stdlib.h"
#include "soil.h" // added by LHE.
#include <math.h>

/********************************************************************/
void soil_water_factor_v2(struct Soil p[])
// output:  dt, fw; LHE. Jan 29. 2013.
{
	double	ft[MAX_LAYERS], fpsisr[MAX_LAYERS];
	double	dtt[MAX_LAYERS];
	double t1, t2;
//	double fw;
	int i;

	t1 = -0.02;
	t2 = 2.0;

	if (p->psim[0] <= 0.000001) // just in case that this function is called before "updatesoilmoisture". LHE
		for (i = 0; i < p->n_layer; i++)
		{
			p->psim[i] = p->psi_sat[i] * pow(p->thetam[i] / p->fei[i], -p->b[i]);
			p->psim[i] = max(p->psi_sat[i], p->psim[i]);  // I see no necessity to use this line unless thetam > fei. LHE May 20, 2015
		}

	for (i = 0; i < p->n_layer; i++)
	{
		if (p->psim[i] > p->psi_min)  // changed 10.0 to psi_min. LHE. Feb. 13, 2013. 
			fpsisr[i] = 1.0 / (1 + pow((p->psim[i] - p->psi_min) / p->psi_min, p->alpha));   //psi_sr in m H2O! This is the old version. LHE. 
		else	
			fpsisr[i] = 1.0;

		if (p->temp_soil_p[i] > 0.0)	// using a previous temperature | soil water factor is calculated before Tm is updated. LHE.
			ft[i] = (1.0 - exp(t1 * pow(p->temp_soil_p[i], t2))); // 1/x 
		else
			ft[i] = 0.0; 

		// ft[i] = min(ft[i], 1.0);  // LHE. 
		//ft[i] = max(0.001, ft[i]); //  changed 0.1 to 0.001. LHE.
		
		fpsisr[i] = fpsisr[i] * ft[i];
	}  

	if (FW_VERSION == 1) 
	{
		for (i = 0; i < p->n_layer; i++) 
			dtt[i] = p->f_root[i] * fpsisr[i]; /* eq. 14 in Ju 2006 */
	}
	else
	{
		for (i = 0; i < p->n_layer; i++) 
			dtt[i] = p->f_root[i];
	}

	//dtt[0] = 0.0;
	double dtt_sum = 0.0;
	for (i = 0; i < p->n_layer; i++) 
		dtt_sum = dtt_sum + dtt[i];/* eq. 14 in JU2006 */

	if (dtt_sum < 0.000001)
		p->f_soilwater = 0.1; // when soil temperatures in all layers are <=0. LHe
	else
	{
		for (i = 0; i < p->n_layer; i++)
		{
			p->dt[i] = dtt[i] / dtt_sum;

	//		if (p->dt[0] < 0.0000001)
	//			printf("%f\n", p->dt[0]);
			if (isnan(p->dt[i]))
				printf("%f\n", p->dt[0]);
		}
		

		//fpsisr[0] = 0;
		double fpsisr_sum = 0;
		for (i = 0; i < p->n_layer; i++)
		{
			fpsisr_sum = fpsisr_sum + fpsisr[i] * p->dt[i]; /* eq. 12, in Chen 2012 GBC; eq 15 in JU */
		}

		p->f_soilwater = max(0.1, fpsisr_sum);
	}


	//printf("fpsisr[0]=%f\n",fpsisr[0]);
	
}

// a very simple fw by LHE July 29, 2014 @home
//double tan_theta = (1.5 - 0.3) / (soil->Layer[2].theta_vfc - soil->Layer[2].theta_vwp);
//fw = 0.3 + (soil->Layer[2].thetam - soil->Layer[2].theta_vwp) * tan_theta;
//fw = max(0.3, fw);

//double tan_theta = (1.0 - 0.3) / (soil->Layer[2].theta_vfc - soil->Layer[2].theta_vwp);
//fw = 0.3 + (soil->Layer[2].thetam - soil->Layer[2].theta_vwp) * tan_theta;
//fw = max(0.3, fw);
//
//fw = min(3.0, fw);
//
//soil->f_soilwater = fw;

