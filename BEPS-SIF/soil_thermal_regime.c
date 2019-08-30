
// Soil thermal regime: update the soil temperature fore each soil layer

// Liming He
// Last update: 09/15/2015

#include "soil.h"
#include <math.h>


void UpdateHeatFlux(struct Soil p[], double Xg_snow, double lambda_snow, double Tsn0, double Tair_annual_mean, double peroid_in_seconds)
{
	int i;

	for (i = 1; i <= p->n_layer; i++)
	{
		if (i < p->n_layer)
			p->G[i] = (p->temp_soil_p[i - 1] - p->temp_soil_p[i]) / (0.5 * p->d_soil[i - 1] / p->lambda[i - 1] + 0.5 * p->d_soil[i] / p->lambda[i]);
		else
			p->G[i] = p->lambda[i-1] * (p->temp_soil_p[i-1] - Tair_annual_mean) / (DEPTH_F + p->d_soil[i-1] * 0.5); // NOTE: sdat->temp: change to annual Tair
		
		if (p->G[i]>200) p->G[i] = 200;              /*juweimin05*/
		if (p->G[i]<-200) p->G[i] = -200;            /*juweimin05*/
	}



	/*juweimin05--------------------------------------------------------------------*/
	for (i = 0; i < p->n_layer; i++) // starting from zero layers. 
	{
		double S = 0;
		p->temp_soil_c[i] = p->temp_soil_p[i] + (p->G[i] - p->G[i + 1] + S) / (p->Cs[i] * p->d_soil[i]) * peroid_in_seconds; // soil temperatures are updated here. LHE
	
		if (p->temp_soil_c[i] > 50.0)  p->temp_soil_c[i] = 50.0; // the two lines are from Ju. 
		if (p->temp_soil_c[i] < -50.0) p->temp_soil_c[i] = -50.0;// I see no necessarity to use this. LHe. May 19, 2015

	}  /* end of ii loop  */

	Update_ice_ratio(p);

	for (i = 0; i < p->n_layer; i++)
		p->temp_soil_p[i] = p->temp_soil_c[i];

}


void Update_Cs(struct Soil p[])
{
	int i;

	for (i = 0; i < p->n_layer; i++)
	{
		// Chen B. (2007) Ecological Modelling 209, 277-300  (equation 18)
		double term1, term2, term3; // LHE. Feb. 12, 2013. 

		term1 = 2.0 * 1.0e+3 * p->density_soil[i] / 2.65;
		term2 = 1.0e+6 * p->thetam[i] * (4.2 * (1 - p->ice_ratio[i]) + 2.09 * p->ice_ratio[i]);
		term3 = 2.5 * 1.0e+6 * p->f_org[i];

		p->Cs[i] = term1 + term2 + term3;
	}
}

// update the frozen status of each soil.

void Update_ice_ratio(struct Soil p[])
{
	int i;
	double Lf0 = 3.34 * 100000; /* latet heat of fusion (liquid: solid) at 0C*/

	for (i = 0; i < p->n_layer; i++)
	{
		// starting to frozen
		if (p->temp_soil_p[i] >= 0.0 && p->temp_soil_c[i] < 0.0 && p->ice_ratio[i] < 1.0)
		{
			// to do: change Gsf to a tmp var.
			double Gsf = (0.0 - p->temp_soil_c[i]) * p->Cs[i] * p->d_soil[i];
			p->ice_ratio[i] = p->ice_ratio[i] + Gsf / Lf0 / 1000.0 / (p->thetam[i] * p->d_soil[i]);
			p->ice_ratio[i] = min(1.0, p->ice_ratio[i]);

			p->temp_soil_c[i] = 0;
		}
		// starting to melt
		else if (p->temp_soil_p[i] <= 0.0 &&  p->temp_soil_c[i] > 0.0 && p->ice_ratio[i] > 0.0)
		{
			double Gsm = (p->temp_soil_c[i] - 0.0) * p->Cs[i] * p->d_soil[i];
			p->ice_ratio[i] = p->ice_ratio[i] - Gsm / Lf0 / 1000.0 / (p->thetam[i] * p->d_soil[i]);
			p->ice_ratio[i] = max(0.0, p->ice_ratio[i]);

			p->temp_soil_c[i] = 0;
		}

		p->ice_ratio[i] = p->ice_ratio[i] * p->thetam_prev[i] / p->thetam[i];
		p->ice_ratio[i] = min(1.0, p->ice_ratio[i]);
	}
}



// Liming HE. Jan 30, 2013.
// Input: thermal_cond, fei, ice_ratio, thetam, kw, ki
// Output: lambda for each layer.
void UpdateSoilThermalConductivity(struct Soil p[])
{
	int i;

	/*to calculate thermal conductivity of each soil layer, advances in water resources 26(2003), 79-93*/
	double ki = 2.1;       /*the thermal conductivity of ice*/
	double kw = 0.61;      /*the thermal conductivity of water*/

	for (i = 0; i < p->n_layer; i++)
	{
		double tmp1, tmp2, tmp3, tmp4;
		tmp1 = pow(p->thermal_cond[i], (1 - p->fei[i])); // dry
		tmp2 = pow(ki, (1.2 * p->thetam[i] * p->ice_ratio[i]) ); // ice.  no source for "1.2"
		tmp3 = pow(kw, p->thetam[i] * (1 - p->ice_ratio[i]) );  // water
		tmp4 = p->thetam[i] / p->fei[i];  // Sr

		p->lambda[i] = (tmp1 * tmp2 * tmp3 - 0.15) * tmp4 + 0.15; // Note: eq. 8. LHE

		p->lambda[i] = max(p->lambda[i], 0.15);     /*juweimin05*/
	}
}
