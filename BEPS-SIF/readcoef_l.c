/*************************************************************************
  Program:     	readcoef.c 
  Description: reading coefficients 
  Written by:   G. Mo  --J. Chen's Group
  Last update:	Dec., 2005 
*****************************************************************************/

#include "beps.h"

void readcoef(lc,stxt,coef)

int short lc;
int stxt;   
double coef[50];
{
  double lignion_leaf,lignion_fr,lignion_wd,clay1,clay_silt1;

  switch (stxt)
   {
	case 1: 
	clay1 = 0.03;
	clay_silt1=0.08;
	break;

	case 2: 
	clay1 = 0.07;
	clay_silt1=0.19;
	break;

	case 3: 
	clay1 = 0.1;
	clay_silt1=0.35;
	break;

	case 4: 
	clay1 = 0.18;
	clay_silt1=0.58;
	break;

	case 5: 
	clay1 = 0.15;
	clay_silt1=0.8;
	break;

	case 6: 
	clay1 = 0.27;
	clay_silt1=0.4;
	break;

	case 7: 
	clay1 = 0.34;
	clay_silt1=0.68;
	break;

	case 8: 
	clay1 = 0.33;
	clay_silt1=0.91;
	break;

	case 9: 
	clay1 = 0.4;
	clay_silt1=0.47;
	break;

	case 10: 
	clay1 = 0.45;
	clay_silt1=0.9;
	break;

	case 11: 
	clay1 = 0.6;
	clay_silt1=0.8;

	case 13: case 14: case 15:
	clay1 = 0.6;
	clay_silt1=0.8; 
	
   }


	switch(lc)
	{
	case 1: case 2: case 3: case 4: case 5: /* conifer */
		coef[0]=0.301013;
		coef[1]=0.148216;
		coef[2]=0.212864;
		coef[3]=0.347907;
		coef[4]=0.02888;
		coef[5]=0.02688;
		coef[6]=0.1925;
		coef[7]=0.5948;
		
		coef[31]=650;
		coef[32]=380;
		coef[33]=100;
		coef[34]=190;
		coef[35]=350.26;
		coef[36]=65.532;
		coef[37]=56.532;
		coef[38]=56.532;
		coef[39]=56.532;
		coef[40]=32.74;
		coef[41]=350.26;
		coef[42]=56.532;
		coef[43]=56.532;
		coef[44]=12;
		coef[45]=10;
		coef[46]=1.3557;
		coef[47]=20;
		
		lignion_leaf=0.302/(0.15+0.018*0.6*coef[36]);
		lignion_fr=0.280/(0.15+0.018*0.6*coef[37]);
		lignion_wd=0.40;
		
		coef[8]=3.9*(exp(-3*lignion_leaf))*((1-lignion_leaf)*0.6+lignion_leaf*0.3);  
		coef[9]=3.9*(exp(-3*lignion_leaf))*(1-lignion_leaf)*0.4;  
		coef[10]=3.9*(exp(-3*lignion_leaf))*lignion_leaf*0.7; 
        
		coef[11]=14.8*0.6;  
		coef[12]=14.8*0.4;  
		
		coef[13]= 4.8*(exp(-3*lignion_fr))*((1-lignion_fr)*0.55+lignion_fr*0.3);  
		coef[14]= 4.8*(exp(-3*lignion_fr))*(1-lignion_fr)*0.45;  
		coef[15]= 4.8*(exp(-3*lignion_fr))*lignion_fr*0.7;  
		
		coef[16]=18.5*0.5;  
		coef[17]=18.5*0.5;  
		
		coef[18]= 2.4*(exp(-3*lignion_wd))*((1-lignion_wd)*0.55+lignion_wd*0.45);  
		coef[19]= 2.4*(exp(-3*lignion_wd))*(1-lignion_wd)*0.45;  
		coef[20]= 2.4*(exp(-3*lignion_wd))*lignion_wd*0.55;    
		
		coef[21]= 7.3*(1-0.75*clay_silt1)*(0.85-0.68*clay_silt1); 
		coef[22]= 7.3*(1-0.75*clay_silt1)*(0.003+0.032*clay1);
		coef[23]= 7.3*(1-0.75*clay_silt1)*(1-(0.003+0.032*clay1)-(0.85-0.68*clay_silt1)-5.0/18.0*(0.01+0.04*(1-clay_silt1))); 
		
		coef[24]=6.0*0.6;
		coef[25]=6.0*0.4;
		
		coef[26]=0.25*0.55;
		coef[27]=0.25*(0.003-0.009*clay1);
		
		if(coef[27]<0.00001)coef[27]=0.00001; 
		if((0.003-0.009*clay1)>0.00001) coef[28]=0.25*(1-0.55-(0.003-0.009*clay1));
		else coef[28]=0.25*(1-0.55-0.00001);
		
		/*    coef[29]=0.0045*0.5;
		coef[30]=0.0045*0.5;*/
		coef[29]=0.007*0.5;
		coef[30]=0.007*0.5;
		break;	

	case 6: case 9: 	/* deciduous, tropic evergreen forest */
		coef[0]=0.422354;
		coef[1]=0.108994;
		coef[2]=0.242626;
		coef[3]=0.226026;
		coef[4]=0.01680;
		coef[5]=0.0248;
		coef[6]=1;
		coef[7]=0.5948;
		coef[31]=650;
		coef[32]=380;
		coef[33]=100;
		coef[34]=190;
		coef[35]=350.26;
		coef[36]=56.532;
		coef[37]=56.532;
		coef[38]=56.532;
		coef[39]=56.532;
		coef[40]=32.74;
		coef[41]=350.26;
		coef[42]=56.532;
		coef[43]=56.532;
		coef[44]=12.0;
		coef[45]=10.0;
		coef[46]=1.3557;
		coef[47]=20.0;
        
		lignion_leaf=0.224/(0.15+0.018*0.6*coef[36]);
		lignion_fr=0.200/(0.15+0.018*0.6*coef[37]);
		lignion_wd=0.30;
		
		coef[8]=3.9*(exp(-3*lignion_leaf))*((1-lignion_leaf)*0.6+lignion_leaf*0.3);  
		coef[9]=3.9*(exp(-3*lignion_leaf))*(1-lignion_leaf)*0.4;  
		coef[10]=3.9*(exp(-3*lignion_leaf))*lignion_leaf*0.7; 
        
		coef[11]=14.8*0.6;  
		coef[12]=14.8*0.4;  
		
		coef[13]= 4.8*(exp(-3*lignion_fr))*((1-lignion_fr)*0.55+lignion_fr*0.3);  
		coef[14]= 4.8*(exp(-3*lignion_fr))*(1-lignion_fr)*0.45;  
		coef[15]= 4.8*(exp(-3*lignion_fr))*lignion_fr*0.7;  
		
		coef[16]=18.5*0.5;  
		coef[17]=18.5*0.5;  
		
		coef[18]= 2.4*(exp(-3*lignion_wd))*((1-lignion_wd)*0.55+lignion_wd*0.45);  
		coef[19]= 2.4*(exp(-3*lignion_wd))*(1-lignion_wd)*0.45;  
		coef[20]= 2.4*(exp(-3*lignion_wd))*lignion_wd*0.55;   
		
		coef[21]= 7.3*(1-0.75*clay_silt1)*(0.85-0.68*clay_silt1); 
		coef[22]= 7.3*(1-0.75*clay_silt1)*(0.003+0.032*clay1);
		coef[23]= 7.3*(1-0.75*clay_silt1)*(1-(0.003+0.032*clay1)-(0.85-0.68*clay_silt1)-5.0/18.0*(0.01+0.04*(1-clay_silt1))); 
		
		coef[24]=6.0*0.6;
		coef[25]=6.0*0.4;
		
		coef[26]=0.25*0.55;
		coef[27]=0.25*(0.003-0.009*clay1);
		
		if(coef[27]<0.00001)coef[27]=0.00001; 
		if((0.003-0.009*clay1)>0.00001) coef[28]=0.25*(1-0.55-(0.003-0.009*clay1));
		else coef[28]=0.25*(1-0.55-0.00001);
		
		coef[29]=0.0045*0.5;
		coef[30]=0.0045*0.5;
		/*    coef[29]=0.007*0.5;
		//    coef[30]=0.007*0.5;*/
		
		break;

	case 10:  	/* mix forest */
		coef[0]=0.4;
		coef[1]=0.1;
		coef[2]=0.3;
		coef[3]=0.2;
		coef[4]=0.03;
		coef[5]=0.03;
		coef[6]=0.5;
		coef[7]=0.5;
		coef[31]=650;
		coef[32]=380;
		coef[33]=100;
		coef[34]=190;
		coef[35]=350.26;
		coef[36]=56.532;
		coef[37]=56.532;
		coef[38]=56.532;
		coef[39]=56.532;
		coef[40]=32.74;
		coef[41]=350.26;
		coef[42]=56.532;
		coef[43]=56.532;
		coef[44]=12.0;
		coef[45]=10.0;
		coef[46]=1.3557;
		coef[47]=20.0;
        
		lignion_leaf=0.224/(0.15+0.018*0.6*coef[36]);
		lignion_fr=0.200/(0.15+0.018*0.6*coef[37]);
		lignion_wd=0.30;
		
		coef[8]=3.9*(exp(-3*lignion_leaf))*((1-lignion_leaf)*0.6+lignion_leaf*0.3);  
		coef[9]=3.9*(exp(-3*lignion_leaf))*(1-lignion_leaf)*0.4;  
		coef[10]=3.9*(exp(-3*lignion_leaf))*lignion_leaf*0.7; 
        
		coef[11]=14.8*0.6;  
		coef[12]=14.8*0.4;  
		
		coef[13]= 4.8*(exp(-3*lignion_fr))*((1-lignion_fr)*0.55+lignion_fr*0.3);  
		coef[14]= 4.8*(exp(-3*lignion_fr))*(1-lignion_fr)*0.45;  
		coef[15]= 4.8*(exp(-3*lignion_fr))*lignion_fr*0.7;  
		
		coef[16]=18.5*0.5;  
		coef[17]=18.5*0.5;  
		
		coef[18]= 2.4*(exp(-3*lignion_wd))*((1-lignion_wd)*0.55+lignion_wd*0.45);  
		coef[19]= 2.4*(exp(-3*lignion_wd))*(1-lignion_wd)*0.45;  
		coef[20]= 2.4*(exp(-3*lignion_wd))*lignion_wd*0.55;   
		
		coef[21]= 7.3*(1-0.75*clay_silt1)*(0.85-0.68*clay_silt1); 
		coef[22]= 7.3*(1-0.75*clay_silt1)*(0.003+0.032*clay1);
		coef[23]= 7.3*(1-0.75*clay_silt1)*(1-(0.003+0.032*clay1)-(0.85-0.68*clay_silt1)-5.0/18.0*(0.01+0.04*(1-clay_silt1))); 
		
		coef[24]=6.0*0.6;
		coef[25]=6.0*0.4;
		
		coef[26]=0.25*0.55;
		coef[27]=0.25*(0.003-0.009*clay1);
		
		if(coef[27]<0.00001)coef[27]=0.00001; 
		if((0.003-0.009*clay1)>0.00001) coef[28]=0.25*(1-0.55-(0.003-0.009*clay1));
		else coef[28]=0.25*(1-0.55-0.00001);
		
		coef[29]=0.0045*0.5;
		coef[30]=0.0045*0.5;
		/*    coef[29]=0.007*0.5;
		//    coef[30]=0.007*0.5;*/
		
		break;

	case 13:  /* shrub  */
		coef[0]=0.189428;
		coef[1]=0.053605;
		coef[2]=0.45;
		coef[3]=0.306967;
		coef[4]=0.025;
		coef[5]=0.04;
		coef[6]=0.8;
		coef[7]=0.75;
	    coef[31]=650;
		coef[32]=380;
		coef[33]=100;
		coef[34]=190;
		coef[35]=370.26;
		coef[36]=63.532;
		coef[37]=63.532;
		coef[38]=63.532;
		coef[39]=63.532;
		coef[40]=32.74;
		coef[41]=370.26;
		coef[42]=63.532;
		coef[43]=63.532;
		coef[44]=12;
		coef[45]=10;
		coef[46]=1.3557;
		coef[47]=20;
		
		lignion_leaf=0.282/(0.15+0.018*0.6*coef[36]);
		lignion_fr=0.24/(0.15+0.018*0.6*coef[37]);
		lignion_wd=0.35;		 
		
		coef[8]=3.9*(exp(-3*lignion_leaf))*((1-lignion_leaf)*0.6+lignion_leaf*0.3);  
        coef[9]=3.9*(exp(-3*lignion_leaf))*(1-lignion_leaf)*0.4;  
        coef[10]=3.9*(exp(-3*lignion_leaf))*lignion_leaf*0.7; 
        
		coef[11]=14.8*0.6;  
        coef[12]=14.8*0.4;  

        coef[13]= 4.8*(exp(-3*lignion_fr))*((1-lignion_fr)*0.55+lignion_fr*0.3);  
        coef[14]= 4.8*(exp(-3*lignion_fr))*(1-lignion_fr)*0.45;  
        coef[15]= 4.8*(exp(-3*lignion_fr))*lignion_fr*0.7;  

        coef[16]=18.5*0.5;  
        coef[17]=18.5*0.5;  

        coef[18]= 2.4*(exp(-3*lignion_wd))*((1-lignion_wd)*0.55+lignion_wd*0.45);  
        coef[19]= 2.4*(exp(-3*lignion_wd))*(1-lignion_wd)*0.45;  
        coef[20]= 2.4*(exp(-3*lignion_wd))*lignion_wd*0.55;    
		 
		coef[21]= 7.3*(1-0.75*clay_silt1)*(0.85-0.68*clay_silt1); 
		coef[22]= 7.3*(1-0.75*clay_silt1)*(0.003+0.032*clay1);
		coef[23]= 7.3*(1-0.75*clay_silt1)*(1-(0.003+0.032*clay1)-(0.85-0.68*clay_silt1)-5.0/18.0*(0.01+0.04*(1-clay_silt1))); 
       	  
        coef[24]=6.0*0.6;
		coef[25]=6.0*0.4;
		
	    coef[26]=0.25*0.55;
		coef[27]=0.25*(0.003-0.009*clay1);

	   if(coef[27]<0.00001)coef[27]=0.00001; 
	     if((0.003-0.009*clay1)>0.00001) coef[28]=0.25*(1-0.55-(0.003-0.009*clay1));
         else coef[28]=0.25*(1-0.55-0.00001);

		coef[29]=0.007*0.5;
		coef[30]=0.007*0.5;
/*		coef[29]=0.0045*0.5;
//		coef[30]=0.0045*0.5;*/
		break;
	case 14: case 15:  /* grassland and cropland */
		coef[0]=0.0;
		coef[1]=0.0;
		coef[2]=0.5;
		coef[3]=0.5;
		coef[4]=0.0001;
		coef[5]=0.0001;
		coef[6]=1.0;
		coef[7]=0.75;
	    coef[31]=650;
		coef[32]=380;
		coef[33]=100;
		coef[34]=190;
		coef[35]=370.26;
		coef[36]=63.532;
		coef[37]=63.532;
		coef[38]=63.532;
		coef[39]=63.532;
		coef[40]=32.74;
		coef[41]=370.26;
		coef[42]=63.532;
		coef[43]=63.532;
		coef[44]=12;
		coef[45]=10;
		coef[46]=1.3557;
		coef[47]=20;
		
		lignion_leaf=0.282/(0.15+0.018*0.6*coef[36]);
		lignion_fr=0.24/(0.15+0.018*0.6*coef[37]);
		lignion_wd=0.35;		 
		
		coef[8]=3.9*(exp(-3*lignion_leaf))*((1-lignion_leaf)*0.6+lignion_leaf*0.3);  
        coef[9]=3.9*(exp(-3*lignion_leaf))*(1-lignion_leaf)*0.4;  
        coef[10]=3.9*(exp(-3*lignion_leaf))*lignion_leaf*0.7; 
        
		coef[11]=14.8*0.6;  
        coef[12]=14.8*0.4;  

        coef[13]= 4.8*(exp(-3*lignion_fr))*((1-lignion_fr)*0.55+lignion_fr*0.3);  
        coef[14]= 4.8*(exp(-3*lignion_fr))*(1-lignion_fr)*0.45;  
        coef[15]= 4.8*(exp(-3*lignion_fr))*lignion_fr*0.7;  

        coef[16]=18.5*0.5;  
        coef[17]=18.5*0.5;  

        coef[18]= 2.4*(exp(-3*lignion_wd))*((1-lignion_wd)*0.55+lignion_wd*0.45);  
        coef[19]= 2.4*(exp(-3*lignion_wd))*(1-lignion_wd)*0.45;  
        coef[20]= 2.4*(exp(-3*lignion_wd))*lignion_wd*0.55;    
		 
		coef[21]= 7.3*(1-0.75*clay_silt1)*(0.85-0.68*clay_silt1); 
		coef[22]= 7.3*(1-0.75*clay_silt1)*(0.003+0.032*clay1);
		coef[23]= 7.3*(1-0.75*clay_silt1)*(1-(0.003+0.032*clay1)-(0.85-0.68*clay_silt1)-5.0/18.0*(0.01+0.04*(1-clay_silt1))); 
       	  
        coef[24]=6.0*0.6;
		coef[25]=6.0*0.4;
		
	    coef[26]=0.25*0.55;
		coef[27]=0.25*(0.003-0.009*clay1);

	   if(coef[27]<0.00001)coef[27]=0.00001; 
	     if((0.003-0.009*clay1)>0.00001) coef[28]=0.25*(1-0.55-(0.003-0.009*clay1));
         else coef[28]=0.25*(1-0.55-0.00001);

		coef[29]=0.007*0.5;
		coef[30]=0.007*0.5;
/*		coef[29]=0.0045*0.5;
//		coef[30]=0.0045*0.5;*/
		break;

	default:  	/* others landcover */
		coef[0]=0.331684;
		coef[1]=0.053605;
		coef[2]=0.307745;
		coef[3]=0.306967;
		coef[4]=0.0278;
		coef[5]=0.0448;
		coef[6]=0.39448;
		coef[7]=0.5948;
		
		coef[31]=650;
		coef[32]=380;
		coef[33]=100;
		coef[34]=190;
		coef[35]=370.26;
		coef[36]=63.532;
		coef[37]=63.532;
		coef[38]=63.532;
		coef[39]=63.532;
		coef[40]=32.74;
		coef[41]=370.26;
		coef[42]=63.532;
		coef[43]=63.532;
		coef[44]=12;
		coef[45]=10;
		coef[46]=1.3557;
		coef[47]=20;
		
		lignion_leaf=0.6*0.224/(0.15+0.018*0.6*coef[36]);
		lignion_fr=0.6*0.200/(0.15+0.018*0.6*coef[37]);
		lignion_wd=0.30;
		
		coef[8]=3.9*(exp(-3*lignion_leaf))*((1-lignion_leaf)*0.6+lignion_leaf*0.3);  
		coef[9]=3.9*(exp(-3*lignion_leaf))*(1-lignion_leaf)*0.4;  
		coef[10]=3.9*(exp(-3*lignion_leaf))*lignion_leaf*0.7; 
        
		coef[11]=14.8*0.6;  
		coef[12]=14.8*0.4;  
		
		coef[13]= 4.8*(exp(-3*lignion_fr))*((1-lignion_fr)*0.55+lignion_fr*0.3);  
		coef[14]= 4.8*(exp(-3*lignion_fr))*(1-lignion_fr)*0.45;  
		coef[15]= 4.8*(exp(-3*lignion_fr))*lignion_fr*0.7;  
		
		coef[16]=18.5*0.5;  
		coef[17]=18.5*0.5;  
		
		coef[18]= 2.4*(exp(-3*lignion_wd))*((1-lignion_wd)*0.55+lignion_wd*0.45);  
		coef[19]= 2.4*(exp(-3*lignion_wd))*(1-lignion_wd)*0.45;  
		coef[20]= 2.4*(exp(-3*lignion_wd))*lignion_wd*0.55;   
		
		coef[21]= 7.3*(1-0.75*clay_silt1)*(0.85-0.68*clay_silt1); 
		coef[22]= 7.3*(1-0.75*clay_silt1)*(0.003+0.032*clay1);
		coef[23]= 7.3*(1-0.75*clay_silt1)*(1-(0.003+0.032*clay1)-(0.85-0.68*clay_silt1)-5.0/18.0*(0.01+0.04*(1-clay_silt1))); 
		
		coef[24]=6.0*0.6;
		coef[25]=6.0*0.4;
		
		coef[26]=0.25*0.55;
		coef[27]=0.25*(0.003-0.009*clay1);
		
		if(coef[27]<0.00001)coef[27]=0.00001; 
		if((0.003-0.009*clay1)>0.00001) coef[28]=0.25*(1-0.55-(0.003-0.009*clay1));
		else coef[28]=0.25*(1-0.55-0.00001);
		
		coef[29]=0.007*0.5;
		coef[30]=0.007*0.5;
		/*    coef[29]=0.0045*0.5;
		//    coef[30]=0.0045*0.5;*/
		
	}

return;
}
