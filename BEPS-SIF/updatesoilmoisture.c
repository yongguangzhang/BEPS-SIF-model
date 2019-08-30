
#include "soil.h"
#include <math.h>

void UpdateSoilMoisture(struct Soil p[], double kstep)
// Last revision: May 20, 2015, by LHe.
// L. HE. Jan 24, 2013.

// Given the current condition, to calcualte soil mositure after a peroid.
// Richards equation. sources: ET and rain

// kstep is defined in beps.h: the total second in this step (the peroid). 
// kkk (outside of the function): step within an hour or half hour measurement. 

{
	double Inf, Inf_max; // infiltration, and Maximum infiltration
	int i;
	 double this_step = 0; //LHE
	 double total_t = 0; //LHE
	 double max_Fb = 0; //LHE
	 
	// assign the current soil temperature to prev variabls. 
	for(i=0; i<=p->n_layer; i++) // save previous thetam. LHE.
		 p->thetam_prev[i] = p->thetam[i];

	for(i=0; i<=p->n_layer; i++)
	{
		 if( p->temp_soil_c[i] > 0.0)  
			  p->f_ice[i] = 1.0; // f_ice should be named as f_water? LHe
		 else if( p->temp_soil_c[i] < -1.0) 
			  p->f_ice[i] = 0.1; // maybe f_ice should be zero? LHE.
		 else p->f_ice[i] = 0.1 + 0.9 * ( p->temp_soil_c[i] + 1.0); 
	}

/*juweimin================================================*/
	 // this part solve the upper boundary  condition(Infiltration). LHE
	 // the maximum Infiltration. Reference? The inf should be changing very fast during a precipaiton because thetam is changing. LHE.
	 Inf_max=p->f_ice[0] * p->Ksat[0] * ( 1 + (p->fei[0] - p->thetam_prev[0]) / p->d_soil[0] * p->psi_sat[0] * p->b[0] / p->fei[0]);

	 Inf = max(p->f_ice[0] * (p->Zp / kstep + p->r_rain_g), 0); // This should be unnessary since tmp >=0; LHE.
	 Inf=min(Inf_max,Inf); 
	 Inf=max(0,Inf);  // This should be unnessary since Inf_max >=0; LHE.

	 p->Zp = (p->Zp / kstep + p->r_rain_g - Inf) * kstep * p->r_drainage; // Ponded water after runoff. This one is related to runoff. LHe. 

/*==============juweimin----------------------------------------*/

	 //// begining of self-adaption step by LHE. Oct 16, 2012
	 while(total_t < kstep) 
	 {
		 for(i=0; i < p->n_layer; i++) 
			 p->km[i] = p->f_ice[i] * p->Ksat[i] * pow((p->thetam[i] / p->fei[i]), (2 * p->b[i] + 3));

		 // soil moisture in the boundaries
		 for(i=0; i < p->n_layer; i++) // 
		 {	
			 if(i < p->n_layer - 1 ) 
				 p->thetab[i] = (p->thetam[i + 1] / p->d_soil[i + 1] + p->thetam[i] / p->d_soil[i]) / (1 / p->d_soil[i] + 1 / p->d_soil[i + 1]);
			 else { // the lowest p->n_layer
				 double d1;
				 d1 = (p->thetam[i] - p->thetab[i-1]) * 2.0 / p->d_soil[i];
				 d1=max(d1,0);
				 p->thetab[i] = p->thetam[i] + d1*p->d_soil[i] / 2.0;
				 p->thetab[i] = min(p->thetab[i], p->fei[i]);
			 }

		 }

		 for(i = 0; i < p->n_layer; i++)
		 {
			 if(i < p->n_layer - 1 )  // the unsaturated hydraulic conductivity at soil lower boundary. 
				 p->Kb[i] = p->f_ice[i] * (p->Ksat[i] * p->d_soil[i] + p->Ksat[i+1] * p->d_soil[i + 1]) / (p->d_soil[i] + p->d_soil[i + 1]) * \
				 pow(p->thetab[i] / p->fei[i], (2 * p->b[i] + 3)); // Note: Kb[0] to Kb[n_layer-1] are not used in the model. LHe.
			 else // when i == p->LAYER - 1
				 p->Kb[i] = 0.5 * p->f_ice[i] * p->Ksat[i] * pow(p->thetab[i] / p->fei[i], (2 * p->b[i] + 3));
		 } 

		 // the unsaturated soil water retention. LHe
		 for(i = 0;i <= p->n_layer - 1;i++) 
		 {
			 p->psim[i] = p->psi_sat[i] * pow(p->thetam[i] / p->fei[i], - p->b[i]);
			 p->psim[i] = max(p->psi_sat[i], p->psim[i] );  // I see no necessity to use this line unless thetam > fei. LHE May 20, 2015

			 // if (p->psim[i] > 300)  p->psim[i] = 300.0;      /*juweimin05   130->300 */
		 }

		 //the unsaturated soil water retention @ boundaries. LHe
		 for(i=0; i < p->n_layer; i++) 
		 {
			 p->psib[i] = p->psi_sat[i] * pow(p->thetab[i] / p->fei[i], - p->b[i]);
			 p->psib[i] = max(p->psi_sat[i], p->psib[i]);
			 // p->psib[i] = min(300,p->psib[i]);
		 } 

		 // the unsaturated hydraulic conductivity of soil p->n_layer @ boundaries
		 for(i=0; i < p->n_layer; i++)
		 {
			 if(i < p-> n_layer-1) 
				 p->KK[i] = (p->km[i] * p->psim[i] + p->km[i + 1] * p->psim[i + 1])  \
							/ (p->psim[i] + p->psim[i + 1]) * (p->b[i] + p->b[i + 1]) / (p->b[i] + p->b[i + 1] + 6);  /*See seller's*/
			 else 	// when i == LAYER-1
				 p->KK[i]  = (p->km[i] * p->psim[i] + p->Kb[i] * p->psib[i])    \
							/(p->psim[i] + p->psib[i]) * p->b[i] / (p->b[i] + 3);
		 }

		 /// // Fb, flow speed. Dancy's law. LHE. 
		 for(i=0;i <= p->n_layer - 1; i++)
		 {
			 if(i < p->n_layer - 1)  
				 {
					 p->r_waterflow[i] = p->KK[i] * (2*(p->psim[i+1] - p->psim[i]) / (p->d_soil[i] + p->d_soil[i+1])+1); /* downwards, positive*/
					 // +1 accounts for gravitational drainage. LHE
				}
			 else
				 //p->r_waterflow[i] = p->km[i] * (0+1); // Seller's 1996. Eq. 37. simplified. 
				 p->r_waterflow[i] = 0; // from Ju. 
		 }

		 
		 // check the r_waterflow further. LHE
		 for(i=0;i<p->n_layer-1;i++)
		 {
				 p->r_waterflow[i] = min((p->fei[i+1] - p->thetam[i+1]) * p->d_soil[i+1] / kstep + p->Ett[i+1], p->r_waterflow[i]); // this line is supposed to confine the flux or flow speed but this is not enough. LHE Oct 16, 2012.
				 if(fabs(p->r_waterflow[i]) > max_Fb) max_Fb = fabs(p->r_waterflow[i]); // find max_Fb for all p->LAYERSs. 
		 }
				 
		 
		 if(max_Fb > 1.0e-5)
			 this_step = 1.0; // determinte the sub-step according to order of Fb empirically .
		 else if(max_Fb > 1.0e-6)
			 this_step = 30.0;
		 else this_step = 360.0;

		 total_t = total_t + this_step;
		 if(total_t > kstep) this_step = this_step - (total_t - kstep);

		 // from there: kstep is replaced by this_step. LHE
		 for(i=0;i < p->n_layer;i++)
		 {
			 if(i==0) 
				 p->thetam[i] = p->thetam[i] + (Inf * this_step - p->r_waterflow[i] * this_step - p->Ett[i] * this_step) / p->d_soil[i]; // kstep->this_step
			 else      
				 p->thetam[i] = p->thetam[i] + (p->r_waterflow[i-1] * this_step - p->r_waterflow[i] * this_step - p->Ett[i] * this_step) / p->d_soil[i];	// kstep->this_step	 

			 /*  thetam[i][kkk]=max(theta_vwp[i]*0.25,thetam[i][kkk]); */
			 p->thetam[i] = max(p->theta_vwp[i], p->thetam[i]); 
			 p->thetam[i] = min(p->fei[i],p->thetam[i]); 
		 }
	
	 } // end of while: the self-adaption step. by LHE.


	 for(i=0;i<p->n_layer;i++)
	 {  // ref?
		 p->ice_ratio[i] = p->ice_ratio[i] * p->thetam_prev[i] / p->thetam[i]; 
		 p->ice_ratio[i] = min(1.0, p->ice_ratio[i]);
	 }

}



// Soil water uptake from a layer. 
// Liming He @ PGB 205B, May 26, 2015
void Soil_Water_Uptake(struct Soil p[], double Trans_o, double Trans_u, double Evap_soil)
{
	int i;
	double rho_w = 1025.0;
	double Source;

	Source = Trans_o + Trans_u;

	// for the top layer
	p->Ett[0] = Source / rho_w * p->dt[0] + Evap_soil / rho_w;

	// for each layer: 
	for (i = 1; i < p->n_layer; i++)
		p->Ett[i] = Source / rho_w * p->dt[i];
}

