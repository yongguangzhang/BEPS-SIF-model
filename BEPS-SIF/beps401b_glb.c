/*************************************************************************
  Program: beps401b_glb.c (main program for global simulation)
	version 4.01b: for simulation of NPP, NEP and respirations
 *************************************************************************/

#include "beps.h"
#include "soil.h"

int main(argc,argv)
int argc;
char *argv[];
{
	int pix,line;		
	int jday,rstep;	
	char argv1[20];
    	short rv,tv,hv,pv,wv,laiv;
    
	int i,j,p,ptr,ptr2;
	int lc,stxt;
	float lon,lat;			
	float **lc0,**lai0,**ci0,**lon0,**lat0; 	
	float **st0,**sw0,**sdp0,**stxt0;
	float **Ccd0,**Cssd0,**Csmd0,**Cfsd0,**Cfmd0,**Csm0,**Cm0,**Cs0,**Cp0,**laiy0,**nppy0;
	
	double ci,lai,laiy,nppy;	

	double tem,hum,st,sw,snowdepth;
	double Ccd[5],Cssd[5],Csmd[5],Cfsd[5],Cfmd[5],Csm[5],Cm[5],Cs[5],Cp[5];

	double parameter[50];
	double coef[100]; 
	double es,esd;
	double CosZs;  

	struct climatedata *meteo;
	struct results *mid_res;
	struct Soil *soilp;

	double var_o[42],var_n[42];
	double v2last[40];
	float outp[10];
	float **ann1,**ann2,**ann3;
	double total[10];

	char outp1fn[80],outp2fn[80],outp3fn[80]; 

	FILE *outp1_p,*outp2_p,*outp3_p;
	FILE *inR_p,*inT_p,*inH_p,*inP_p,*wind_p,*laif_p;

/** Reading control file **/
	readconf(argv[1]);

	printf("running with %s  ...\n",argv[1]);
	
	meteo=(struct climatedata*)malloc(366*sizeof(struct climatedata));
	mid_res=(struct results*)malloc(366*sizeof(struct results));
	soilp = (struct Soil*)malloc(sizeof(struct Soil));
		

	lc0 = (float **) malloc((nlines*2)*sizeof(float*));
		for (j=0;j<=nlines+1;j++)
          lc0[j] = (float *) malloc((npixels+1)*sizeof(float));
    	ci0 = (float **) malloc((nlines*2)*sizeof(float*));
           for (j=0;j<=nlines+1;j++)
          ci0[j] = (float *) malloc((npixels+1)*sizeof(float));
    	lon0 = (float **) malloc(nlines*2*sizeof(float*));
		for (j=0;j<=nlines;j++)
          lon0[j] = (float *) malloc(npixels*sizeof(float));
    	lat0 = (float **) malloc(nlines*2*sizeof(float*));
		for (j=0;j<=nlines;j++)
          lat0[j] = (float *) malloc(npixels*sizeof(float));
	stxt0 = (float **) malloc((nlines*2)*sizeof(float*));
		for (j=0;j<=nlines+1;j++)
          stxt0[j] = (float *) malloc((npixels+1)*sizeof(float));
    	st0 = (float **) malloc((nlines*2)*sizeof(float*));
		for (j=0;j<=nlines+1;j++)
          st0[j] = (float *) malloc((npixels+1)*sizeof(float));
    	sw0 = (float **) malloc((nlines*2)*sizeof(float*));
		for (j=0;j<=nlines+1;j++)
          sw0[j] = (float *) malloc((npixels+1)*sizeof(float));
    	sdp0 = (float **) malloc((nlines*2)*sizeof(float*));
		for (j=0;j<=nlines+1;j++)
          sdp0[j] = (float *) malloc((npixels+1)*sizeof(float));

    	nppy0 = (float **) malloc(361*sizeof(float));
		for (j=0;j<=181;j++)
          nppy0[j] = (float *) malloc(361*sizeof(float));
    	laiy0 = (float **) malloc(361*sizeof(float));
           for (j=0;j<=181;j++)
          laiy0[j] = (float *) malloc(361*sizeof(float));

    	Ccd0 = (float **) malloc(361*sizeof(float));
	   	for (j=0;j<=181;j++)
          Ccd0[j] = (float *) malloc(361*sizeof(float));
    	Cssd0 = (float **) malloc(361*sizeof(float));
		for (j=0;j<=181;j++)
          Cssd0[j] = (float *) malloc(361*sizeof(float));
   	Csmd0 = (float **) malloc(361*sizeof(float));
		for (j=0;j<=181;j++)
          Csmd0[j] = (float *) malloc(361*sizeof(float));
    	Cfsd0 = (float **) malloc(361*sizeof(float));
		for (j=0;j<=181;j++)
          Cfsd0[j] = (float *) malloc(361*sizeof(float));
    	Cfmd0 = (float **) malloc(361*sizeof(float));
		for (j=0;j<=181;j++)
          Cfmd0[j] = (float *) malloc(361*sizeof(float));
    	Csm0 = (float **) malloc(361*sizeof(float));
		for (j=0;j<=181;j++)
          Csm0[j] = (float *) malloc(361*sizeof(float));
    	Cm0 = (float **) malloc(361*sizeof(float));
		for (j=0;j<=181;j++)
          Cm0[j] = (float *) malloc(361*sizeof(float));
    	Cs0 = (float **) malloc(361*sizeof(float));
		for (j=0;j<=181;j++)
          Cs0[j] = (float *) malloc(361*sizeof(float));
    	Cp0 = (float **) malloc(361*sizeof(float));
		for (j=0;j<=181;j++)
          Cp0[j] = (float *) malloc(361*sizeof(float));
    	lai0 = (float **) malloc(760*sizeof(float));
		for (j=0;j<=375;j++)
          lai0[j] = (float *) malloc(73000*sizeof(float));

	ann1 = (float **) malloc((nlines*2)*sizeof(float*));
		for (j=0;j<=nlines+1;j++)
          ann1[j] = (float *) malloc((npixels+1)*sizeof(float));
	ann2 = (float **) malloc((nlines*2)*sizeof(float*));
		for (j=0;j<=nlines+1;j++)
          ann2[j] = (float *) malloc((npixels+1)*sizeof(float));
	ann3 = (float **) malloc((nlines*2)*sizeof(float*));
		for (j=0;j<=nlines+1;j++)
          ann3[j] = (float *) malloc((npixels+1)*sizeof(float));




/*****   Reading input data set 1  *****/
	readinput1(lc0,ci0,stxt0,st0,sw0,sdp0);

/*****   Reading long. lat. values  *****/
	readlonlat(lon0,lat0);

/*****   Reading input dataset2 for calulation of respirations  *****/
	readinput2(laiy0,nppy0,Ccd0,Cssd0,Csmd0,Cfsd0,Cfmd0,Csm0,Cm0,Cs0,Cp0);

/*****   Read lai data  */
	readlai_all(lai0);

/* open climate data files */
	if ((inR_p=fopen(r_fn, "rb")) == NULL) 
	{
       	printf("\n Unable to open r file,  exitting ...\n",r_fn);
       	exit(0);
    	}
	if ((inT_p=fopen(t_fn, "rb")) == NULL) 
	{
       	printf("\n Unable to open file <%s>,  exitting ...\n",t_fn);
       	exit(0);
    	}
	if ((inH_p=fopen(h_fn, "rb")) == NULL) 
	{
       	printf("\n Unable to open file <%s>, exitting ...\n",h_fn);
       	exit(0);
    	}
	if ((inP_p=fopen(p_fn, "rb")) == NULL) 
	{
       	printf("\n Unable to open file <%s>,  exitting ...\n",p_fn);
       	exit(0);
    	}
	if ((wind_p=fopen(wd_fn, "rb")) == NULL) 
	{
       	printf("\n Unable to open file <%s>,  exitting ...\n",wd_fn);
       	exit(0);
    	}


          /*    Open hourly output file */
	/*sprintf(outp1fn,"%s_%dhr.bin",fp4outp1,yr);
    printf("output file name: %s",outp1fn);
    sprintf(outp2fn,"%s_%dhr.bin",fp4outp2,yr);
    sprintf(outp3fn,"%s_%dhr.bin",fp4outp3,yr);

        if ((outp1_p=fopen(outp1fn, "wb")) == NULL)
        {   printf("\n Unable to open output1, exitting...\n");
                exit(0);
        }
        if ((outp2_p=fopen(outp2fn, "wb")) == NULL)
        {   printf("\n Unable to open output2, exitting ...\n");
                exit(0);
        }
        if ((outp3_p=fopen(outp3fn, "wb")) == NULL)
        {   printf("\n Unable to open output3, exitting ...\n");
                exit(0);
        }
    */



	for (line=0; line <nlines; line++)
	{ 
	  printf("line = %d \n", line);
	  for (pix=0; pix<npixels; pix++)
	  {
		lon=lon0[line][pix];
		lat=lat0[line][pix];

		lc = (int)lc0[line][pix];
		ci=ci0[line][pix];

		stxt=(int)stxt0[line][pix];
		st=st0[line][pix];
		sw=sw0[line][pix];
		snowdepth=sdp0[line][pix];

		laiy=laiy0[line][pix];
		nppy=nppy0[line][pix];

		Ccd[0]=Ccd0[line][pix]*1000;   // unit: kg to g
		Cssd[0]=Cssd0[line][pix]*1000;
		Csmd[0]=Csmd0[line][pix]*1000;
		Cfsd[0]=Cfsd0[line][pix]*1000;
		Cfmd[0]=Cfmd0[line][pix]*1000;
		Csm[0]=Csm0[line][pix]*1000;
		Cm[0]=Cm0[line][pix]*1000;
		Cs[0]=Cs0[line][pix]*1000;
		Cp[0]=Cp0[line][pix]*1000;
		if (lc == 9)    //evgn
	    {
	        Ccd[0]=Ccd[0]*1.3;   
		    Cssd[0]=Cssd[0]*1.3;
		    Csmd[0]=Csmd[0]*1.3;
		    Cfsd[0]=Cfsd[0]*1.3;
		    Cfmd[0]=Cfmd[0]*1.3;
		    Csm[0]=Csm[0]*1.3;
		    Cm[0]=Cm[0]*1.3;
		    Cs[0]=Cs[0]*1.3;
		    Cp[0]=Cp[0]*1.3;
	    }
	    if (lc == 13)     //shrub
	    {
	        Ccd[0]=Ccd[0]*1.6;   
		    Cssd[0]=Cssd[0]*1.6;
		    Csmd[0]=Csmd[0]*1.6;
		    Cfsd[0]=Cfsd[0]*1.6;
		    Cfmd[0]=Cfmd[0]*1.6;
		    Csm[0]=Csm[0]*1.6;
		    Cm[0]=Cm[0]*1.6;
		    Cs[0]=Cs[0]*1.6;
		    Cp[0]=Cp[0]*1.6;
	    }
	    if (lc == 6)     //deci
	    {
	        Ccd[0]=Ccd[0]*1.0;   
		    Cssd[0]=Cssd[0]*1.0;
		    Csmd[0]=Csmd[0]*1.0;
		    Cfsd[0]=Cfsd[0]*1.0;
		    Cfmd[0]=Cfmd[0]*1.0;
		    Csm[0]=Csm[0]*1.0;
		    Cm[0]=Cm[0]*1.0;
		    Cs[0]=Cs[0]*1.0;
		    Cp[0]=Cp[0]*1.0;
	    }

	    if (lc == 10)  //mix forest
	    {
	        Ccd[0]=Ccd[0]*0.93;   
		    Cssd[0]=Cssd[0]*0.93;
		    Csmd[0]=Csmd[0]*0.93;
		    Cfsd[0]=Cfsd[0]*0.93;
		    Cfmd[0]=Cfmd[0]*0.93;
		    Csm[0]=Csm[0]*0.93;
		    Cm[0]=Cm[0]*0.93;
		    Cs[0]=Cs[0]*0.93;
		    Cp[0]=Cp[0]*0.93;
	    }
	    if (lc == 15)   //crop land
	    {
	        Ccd[0]=Ccd[0]*1.15;   
		    Cssd[0]=Cssd[0]*1.15;
		    Csmd[0]=Csmd[0]*1.15;
		    Cfsd[0]=Cfsd[0]*1.15;
		    Cfmd[0]=Cfmd[0]*1.15;
		    Csm[0]=Csm[0]*1.15;
		    Cm[0]=Cm[0]*1.15;
		    Cs[0]=Cs[0]*1.15;
		    Cp[0]=Cp[0]*1.15;
	    }

	    if (lc == 14)   //grassland
	    {
	        Ccd[0]=Ccd[0]*1.5;   
		    Cssd[0]=Cssd[0]*1.5;
		    Csmd[0]=Csmd[0]*1.5;
		    Cfsd[0]=Cfsd[0]*1.5;
		    Cfmd[0]=Cfmd[0]*1.5;
		    Csm[0]=Csm[0]*1.5;
		    Cm[0]=Cm[0]*1.5;
		    Cs[0]=Cs[0]*1.5;
		    Cp[0]=Cp[0]*1.5;
	    }
	    if (lc == 1)   //coni
	    {
	        Ccd[0]=Ccd[0]*0.8;   
		    Cssd[0]=Cssd[0]*0.8;
		    Csmd[0]=Csmd[0]*0.8;
		    Cfsd[0]=Cfsd[0]*0.8;
		    Cfmd[0]=Cfmd[0]*0.8;
		    Csm[0]=Csm[0]*0.8;
		    Cm[0]=Cm[0]*0.8;
		    Cs[0]=Cs[0]*0.8;
		    Cp[0]=Cp[0]*0.8;
	    }

	    if (lc == 2)   //cond
	    {
	        Ccd[0]=Ccd[0]*0.8;   
		    Cssd[0]=Cssd[0]*0.8;
		    Csmd[0]=Csmd[0]*0.8;
		    Cfsd[0]=Cfsd[0]*0.8;
		    Cfmd[0]=Cfmd[0]*0.8;
		    Csm[0]=Csm[0]*0.8;
		    Cm[0]=Cm[0]*0.8;
		    Cs[0]=Cs[0]*0.8;
		    Cp[0]=Cp[0]*0.8;
	    }








		/* Read parameters according to land cover type */
		readparam(lc,parameter);

		readcoef(lc,stxt,coef);  // for soil repiration and NEP calculation 

		total[1]=0;
		total[2]=0; 
		total[3]=0; 

	    /* ---  day loop  --- */
  	    for (jday=bgn_day; jday<=end_day; jday++) 
  	    {
	     
	    	/* re-calc LAI & renew clump index */
		p=line*npixels+pix;
	    	lai=lai0[jday][p]*parameter[2]/ci; 

		for (rstep=0;rstep<24;rstep++) 
		{
		  fread(&rv,sizeof(short),1,inR_p);
		  fread(&tv,sizeof(short),1,inT_p);
		  fread(&hv,sizeof(short),1,inH_p);
		  fread(&pv,sizeof(short),1,inP_p);
		  fread(&wv,sizeof(short),1,wind_p);

		  meteo->Srad = rv*0.1;
		  meteo->temp = tv*0.1;
		  meteo->rain = pv*0.0001;
		  meteo->wind = wv*0.1;
		  meteo->LR = -200.0;   //  -200.0 means no measured longwave radiation, the value will be calculated later 
	
		  tem = tv*0.1;
		  hum = hv*0.1;

		/* vapour pressure in mbar */ 
		  //es =0.46*hum*(tem+273.16)/100;	
		  //esd  =  6.1078 * exp((17.269*tem)/(237.3 + tem));
		/* Calculating Relative humidity */
		  //if (es/esd>=1)  meteo->rh = 100;
		  //else  meteo->rh = 100*es/esd; 

		  meteo->rh = hum;  // when reading in relative humidity

	     if (lc<=0 || sw<=0 || stxt<=0 || ci<=0.5)  // pix without enough data
	     {	
		  outp[1] = 0;
		  outp[2] = 0; 
		  outp[3] = 0;  
	     }
	     else   // starting simulation
	     { 
		if(jday==bgn_day && rstep==0)  // for 1st time step, to initialize var.
		{
		  Init_Soil_Parameters(lc, stxt, parameter[27], soilp); 
		  soilp->r_drainage = parameter[26];
		  Init_Soil_Status(soilp, st, tem, sw, snowdepth); 

		  for (i=0;i<=40;i++)   var_o[i] = 0;
		  for (i=3;i<=8;i++)   var_o[i] = tem;

		  for(i=9;i<=14;i++) var_o[i] = soilp->temp_soil_p[i-9];
		  for(i=21;i<=26;i++) var_o[i] = soilp->thetam_prev[i-21];
		  for(i=27;i<=32;i++) var_o[i] = soilp->ice_ratio[i-27];
		}
		else   //  for other time step
		  for (i=0;i<=40;i++)   var_o[i] = v2last[i];
		
		/* ---  Calculating cos_solar zenith angle Z  --- */
		  s_coszs(jday,rstep,lat,lon,&CosZs); 

		/* --- Running simulation modules   --- */

		  inter_prg(jday,rstep,lai,ci,parameter,meteo,CosZs,var_o,var_n,soilp,mid_res);

		  for (i=0;i<=40;i++)  v2last[i]=var_n[i];


		/*  ---  plant respiration /NPP module  ---  */
		
		plantresp(lc,mid_res,laiy,lai,tem,soilp->temp_soil_c[1],CosZs); 
			   
		/*  ---  soil respiration module  ---  */

		soilresp(Ccd,Cssd,Csmd,Cfsd,Cfmd,Csm,Cm,Cs,Cp,nppy,coef,stxt,soilp,mid_res); 


		  outp[1] = (float)(mid_res->GPP);
		  //outp[2]=(float)(mid_res->gpp_o_sunlit+mid_res->gpp_u_sunlit)*12*step*0.000001; 
		  //outp[3]=(float)(mid_res->gpp_o_shaded+mid_res->gpp_u_shaded)*12*step*0.000001;
		  outp[2] = (float)(mid_res->NPP);
		  outp[3] = (float)(mid_res->NEP);   

	      }   //  end of simulation  if lc<0 else  


		/* ---  Writing outputs of this time step --- */

           	  //fwrite(&outp[1],sizeof(float),1,outp1_p);
           	  //fwrite(&outp[2],sizeof(float),1,outp2_p);
          	  //fwrite(&outp[3],sizeof(float),1,outp3_p);
	
		  total[1]=total[1]+outp[1];
		  total[2]=total[2]+outp[2]; 
		  total[3]=total[3]+outp[3];   

		}   /*  end of rstep - hourly cycle  */
	    }  /* end of for day loop */

	  ann1[line][pix]=total[1];
	  ann2[line][pix]=total[2];
	  ann3[line][pix]=total[3];

	  }   /* end of line/pix loop */
	}

	//fclose(outp1_p);
	//fclose(outp2_p);
	//fclose(outp3_p);

	
	/* ---  Writing annual outputs to files (optional) --- */

	sprintf(outp1fn,"%s_%dann.bin",fp4outp1,yr);
    sprintf(outp2fn,"%s_%dann.bin",fp4outp2,yr);
    sprintf(outp3fn,"%s_%dann.bin",fp4outp3,yr);

        if ((outp1_p=fopen(outp1fn, "wb")) == NULL)
        {   printf("\n Unable to open ann. output1, exitting...\n");
                exit(0);
        }
        if ((outp2_p=fopen(outp2fn, "wb")) == NULL)
        {   printf("\n Unable to open ann. output2, exitting ...\n");
                exit(0);
        }
        if ((outp3_p=fopen(outp3fn, "wb")) == NULL)
        {   printf("\n Unable to open ann. output3, exitting ...\n");
                exit(0);
        }
	

	for (line=0; line<nlines; line++) 
	  for (pix=0; pix<npixels; pix++)
	  {
           fwrite(&ann1[line][pix],sizeof(float),1,outp1_p);
           fwrite(&ann2[line][pix],sizeof(float),1,outp2_p);
           fwrite(&ann3[line][pix],sizeof(float),1,outp3_p);
	  }

	fclose(outp1_p);
	fclose(outp2_p);
	fclose(outp3_p);
       

  return 1;
}	/* end of main */
