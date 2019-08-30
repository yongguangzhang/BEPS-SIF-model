/*************************************************************************
Program:     plant_respir.c
--------
Description:	Estimate plant respiration.

  Written by:   J. Liu. and W.Ju,  Modified by G. Mo
  Last update:  May 2015
  
*****************************************************************************/

#include "beps.h" 

void plantresp(LC,mid_res,lai_yr,lai,temp_air,temp_soil,CosZs) 

int LC;
struct results mid_res[];
double lai_yr;
double temp_air;
double lai,temp_soil;
double CosZs;

{
	double temp_opt25=25.0;
	double biomass,biomass_leaf_o,biomass_stem_o,biomass_root_o,biomass_leaf_u,biomass_stem_u,biomass_root_u;
	double respir_croot_o,respir_root_o,respir_stem_o,respir_leaf_o;
	double respir_croot_u,respir_root_u,respir_stem_u,respir_leaf_u;
	double q10;	
	double exponent;
	double lai_u,lai_max_o,lai_max_u;   
	double ra;
	double coef_leaf_respir,coef_stem_respir,coef_root_respir,coef_fineroot_respir;
	double gpp_o,gpp_u,gpp_r,rg,ratio_froot;
	
	
	if (LC==25 || LC==40) lai_u=0.01;
	else 
		lai_u=1.18*exp(-0.99*lai);
	
	if (lai_u>lai) lai_u=0.01;
	
	if(LC==6)
		ra=0.6;
	else
		ra=1.0;
	
	q10=3.22-0.046*temp_air;
	
	
	switch(LC)
	{
	case 1: case 2: case 3: case 4: case 5: /* conifer */
		/*  calculating aboveground biomass based on LAI  J. Liu 2002 */
		biomass=0.9097*lai_yr+0.125*lai_yr*lai_yr;
		biomass_leaf_o=0.05*biomass;	/* leaf C of overstory */
		biomass_stem_o=0.95*biomass;	/* stem C of overstory */
		biomass_root_o=0.454*biomass;
		/*biomass_root_o=0.232*biomass; // root C of overstoryKurz 1996 */
		biomass_leaf_u=0.3*biomass_leaf_o;	/* leaf C of understory */
		biomass_stem_u=0.02*biomass_stem_o;  	/* stem C of understory */
		biomass_root_u=0.05*biomass_root_o;	/* root C of understory */

		coef_leaf_respir=0.0015/RTIMES;  /*  leaf resp co.  kg C-1 d-1 kg-1    */
        	coef_stem_respir=0.0020/RTIMES;  /*  stem resp co.   kg C-1 d-1 kg-1      */
        	coef_root_respir=0.0020/RTIMES;  /*  root resp co.   kg C-1 d-1 kg-1   */
		coef_fineroot_respir=0.003/RTIMES;    /*  fine root resp co.   kg C-1 d-1 kg-1	*/

		lai_max_o=4.5;		/* LAI_max_overstory                         3.3*/
		lai_max_u=2.4;		/* LAI_max_understroy                        2.4*/

		
		break;	
		
	case 6: case 9: 	/* deciduous forest, tropic evergreen */
		/*  calculating aboveground biomass based on LAI  J. Liu 2002 */
		biomass=1.545*lai_yr+0.183*lai_yr*lai_yr;
		biomass_leaf_o=0.04*biomass;	/* leaf C of overstory */
		biomass_stem_o=0.96*biomass;	/* stem C of overstory */
		biomass_root_o=1.432*pow(biomass,0.639);	/* root C of overstory  Kurz 1996 */
		biomass_leaf_u=0.3*biomass_leaf_o;	/* leaf C of understory */
		biomass_stem_u=0.01*biomass_stem_o;  	/* stem C of understory */
		biomass_root_u=0.01*biomass_root_o;	/* root C of understory */

		coef_leaf_respir=0.015/RTIMES;  /*  leaf resp co.  kg C-1 d-1 kg-1    */
        	coef_stem_respir=0.0035/RTIMES;  /*  stem resp co.   kg C-1 d-1 kg-1      */
        	coef_root_respir=0.0025/RTIMES;  /*  root resp co.   kg C-1 d-1 kg-1   */
		coef_fineroot_respir=0.003/RTIMES;    /*  fine root resp co.   kg C-1 d-1 kg-1	*/

		lai_max_o=4.5;		/* LAI_max_overstory                         3.3*/
		lai_max_u=2.4;		/* LAI_max_understroy                        2.4*/

		break;
		
	case 13:  /* shrub  */
		biomass=1.545*lai_yr+0.183*lai_yr*lai_yr;
		biomass_leaf_o=0.1*biomass;	/* leaf C of overstory */
		biomass_stem_o=0.90*biomass;	/* stem C of overstory */
		biomass_root_o=1.432*pow(biomass,0.639);	/* root C of overstory  Kurz 1996 */
		biomass_leaf_u=0.3*biomass_leaf_o;	/* leaf C of understory */
		biomass_stem_u=0.01*biomass_stem_o;  	/* stem C of understory */
		biomass_root_u=0.01*biomass_root_o;	/* root C of understory */

		coef_leaf_respir=0.001/RTIMES;  /*  leaf resp co.  kg C-1 d-1 kg-1    */
        	coef_stem_respir=0.002/RTIMES;  /*  stem resp co.   kg C-1 d-1 kg-1      */
        	coef_root_respir=0.0015/RTIMES;  /*  root resp co.   kg C-1 d-1 kg-1   */
		coef_fineroot_respir=0.003/RTIMES;    /*  fine root resp co.   kg C-1 d-1 kg-1	*/

		lai_max_o=3.3;  /*   LAI_max_overstory   */
		lai_max_u=0.01; /*  LAI_max_understroy   */

		break;
		
	case 25: case 40:
		biomass_leaf_o=0.05*lai_yr;	/* leaf C = lai/20  from W.Ju 05y11*/
		biomass_stem_o=0.0;	/* stem C */
		biomass_root_o=0.061*lai_yr;	/* root C = lai/20*0.55/0.45  from W.Ju 05y11*/
		biomass_leaf_u=0.0;	
		biomass_stem_u=0.0;  	
		biomass_root_u=0.0;

		coef_leaf_respir=0.001/RTIMES;  /* leaf resp co.  kg C-1 d-1 kg-1    */
        	coef_stem_respir=0.002/RTIMES;  /*  stem resp co.   kg C-1 d-1 kg-1      */
        	coef_root_respir=0.0015/RTIMES;  /* root resp co.   kg C-1 d-1 kg-1   */
		coef_fineroot_respir=0.003/RTIMES;    /*  fine root resp co.   kg C-1 d-1 kg-1	*/

		lai_max_o=3.3;  /*   LAI_max_overstory     */
		lai_max_u=0.01; /*  LAI_max_understroy     */

	}   /* end of switch  */
	
	/* ---- calculation for overstory  ---- */
	
	/* stem maintenance respiration */
	exponent=(temp_air-temp_opt25)/10.0;    
	respir_stem_o = (biomass_stem_o*0.35/(biomass_stem_o+0.35))*coef_stem_respir*pow(q10,exponent)*ra;
	respir_stem_o = max(respir_stem_o, 0.0);
	
	/* root maintenance respiration     changed in Nov.2005 by W. Ju */
	/*	exponent=(temp_soil-temp_opt25)/10.0;
	respir_root_o = biomass_root_o*coef_root_respir*pow(q10,exponent)*ra;
	respir_root_o = max  (respir_root_o, 0.0);   */   
	exponent=(temp_soil-temp_opt25)/10.0;
	
	if (LC==25 || LC==40) 
		respir_root_o = biomass_root_o*coef_root_respir*pow(q10,exponent)*ra;
	else
	{
		ratio_froot=exp(1.007)*pow(biomass_root_o,-(0.841));
		ratio_froot=min(0.9, ratio_froot);
		
		respir_croot_o=0.05*biomass_root_o*(1-ratio_froot)*coef_root_respir*pow(q10,exponent);		/*coarse root */
		respir_root_o=respir_croot_o+0.05*biomass_root_o*ratio_froot*coef_fineroot_respir*pow(q10,exponent);	/* coarse + fine root */
	}
	respir_root_o = max  (respir_root_o, 0.0);  
	
	
	/* leaf day/night time maintenance respiration */
	if (CosZs>0.01) respir_leaf_o=0;
	else 
	{  
		exponent=(temp_air-temp_opt25)/10.0;
		respir_leaf_o =lai/lai_max_o*biomass_leaf_o*coef_leaf_respir*pow(q10,exponent)*ra;
	} 
	respir_leaf_o =max( respir_leaf_o, 0.0);
	
	/*     changed in Nov.2005 by W. Ju */
	//  *npp_o = 0.8*gpp_o-respir_leaf_o-respir_stem_o-respir_root_o; 
        
	gpp_o = (mid_res->gpp_o_sunlit + mid_res->gpp_o_shaded)*12*step*0.000001;            
	gpp_r = gpp_o - (respir_leaf_o+respir_stem_o+respir_root_o)*1000;                  	
	if (gpp_r <= 0)  rg = 0;
	else rg = 0.35*gpp_r;
	
	mid_res->npp_o = (gpp_r - rg);   //NPP_overstory;  gC/m2
	
	/* ---- calculation for understory  ---- */
	
	/* stem maintenance respiration */
	exponent=(temp_air-temp_opt25)/10.0;
	respir_stem_u =(biomass_stem_u*0.35/(biomass_stem_u+0.35))*coef_stem_respir*pow(q10,exponent)*ra;
	respir_stem_u = max(respir_stem_u, 0.0);
	
	/* root maintenance respiration      changed in Nov.2005 by W. Ju */
    exponent=(temp_soil-temp_opt25)/10.0;
	if (LC==25 || LC==40) 
		respir_root_u = biomass_root_u*coef_root_respir*pow(q10,exponent)*ra;
	
	else
	{
		ratio_froot=exp(1.007)*pow(biomass_root_u,-(0.841));
		ratio_froot=min(0.9, ratio_froot);
		
		respir_croot_u=0.05*biomass_root_u*(1-ratio_froot)*coef_root_respir*pow(q10,exponent);		/*coarse root */
		respir_root_u=respir_croot_u+0.05*biomass_root_u*ratio_froot*coef_fineroot_respir*pow(q10,exponent);	/* coarse + fine root */
	}
	
	respir_root_u = max(respir_root_u, 0.0); 
	

	if (CosZs>0.01) respir_leaf_u=0;
	else 
	{
		exponent=(temp_air-temp_opt25)/10.0;
		respir_leaf_u =lai_u/lai_max_u*biomass_leaf_u*coef_leaf_respir*pow(q10,exponent)*ra*0.5;
	}
	respir_leaf_u =max(respir_leaf_u, 0.0);
	
	/*     changed in Nov.2005 by W. Ju */
	//	*npp_u = 0.8*gpp_u-respir_leaf_u-respir_stem_u-respir_root_u; 

	gpp_u = (mid_res->gpp_u_sunlit + mid_res->gpp_u_shaded)*12*step*0.000001;                             	
	gpp_r = gpp_u - (respir_leaf_u+respir_stem_u+respir_root_u)*1000;                  	
	if (gpp_r <= 0)  rg = 0;
	else rg = 0.35*gpp_r;
	
	mid_res->npp_u = (gpp_r - rg); //understory NPP  gC/m2
	

}



