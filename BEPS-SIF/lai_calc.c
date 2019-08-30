#include "beps.h"


//Bin Chen: 30/05/2015

void lai2(stem_o,stem_u,LC,CosZs,GMU,lai_o,clumping,lai_u,lai_o_sunlit,lai_o_shaded,lai_u_sunlit,lai_u_shaded,PAI_o_sunlit,PAI_o_shaded,PAI_u_sunlit,PAI_u_shaded)

double stem_o; // overstory woody area
double stem_u; // understory woody area
int LC; // landcover type
double CosZs; // cosine solar zenith angle
double lai_o; // overstory lai
double GMU;
double clumping; // clumping index
double lai_u; // understory lai
double *lai_o_sunlit; // overstory sunlit lai
double *lai_o_shaded; // overstory shaded lai
double *lai_u_sunlit ; // understory sunlit lai
double *lai_u_shaded; // understory shaded lai

double *PAI_o_sunlit ; // overstory sunlit lai
double *PAI_o_shaded ; // overstory shaded lai
double *PAI_u_sunlit ; // understory sunlit lai
double *PAI_u_shaded ; // understory shaded lai

{
	 if(CosZs>0)
		 *PAI_o_sunlit =2*CosZs*(1-exp(-GMU*clumping*(lai_o+stem_o)/CosZs));
	 else
		 *PAI_o_sunlit =0;

	 *PAI_o_shaded =(lai_o+stem_o)-*PAI_o_sunlit ;
	 
	 if(CosZs>0)
		 *PAI_u_sunlit =2*CosZs*(1-exp(-GMU*clumping*(lai_o+stem_o+lai_u+stem_u)/CosZs))-*PAI_o_sunlit ;
	 else
		 *PAI_u_sunlit =0;
	 
	 *PAI_u_shaded =(lai_u+stem_u)-*PAI_u_sunlit ;
	 
	 if(CosZs>0)
		 *lai_o_sunlit=2*CosZs*(1-exp(-GMU*clumping*lai_o/CosZs));
	 else
		 *lai_o_sunlit=0;
	 
	 *lai_o_shaded=max(0,lai_o-*PAI_o_sunlit );
	 
	 if(CosZs>0)
		 *lai_u_sunlit =2*CosZs*(1-exp(-GMU*clumping*(lai_o+lai_u)/CosZs))-*lai_o_sunlit;
	 else
		 *lai_u_sunlit =0;
	 
	 *lai_u_shaded=max(0,lai_u-*PAI_u_sunlit );	 


}
