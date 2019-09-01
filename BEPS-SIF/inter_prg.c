/*************************************************************************
Program:      inter_prg.c
--------
Description:  a inter-program between main program and modules 
Last update:	July, 2015

Description:  Solar-induced chlorophyll fluorescence at far-red band module
              was coupled.  Bo Qiu
Last update:	Aug, 2019.
*****************************************************************************/

#include "beps.h"
#include "soil.h"

void   inter_prg(jday,rstep,lai,clumping,parameter,meteo,CosZs,var_o,var_n,soilp,mid_res) 

int jday,rstep;
double lai,clumping;
double parameter[];
struct climatedata meteo[];
double CosZs;
double var_o[],var_n[];
struct Soil soilp[];
struct results mid_res[];
{
	int num,kkk,i;
	int landcover;
	double lai_o,lai_u;
	double stem_o,stem_u;

	double d_soil[layer+1];
	double Zsp;                /*the depth of snow on the surface */
	double Zp,Zp1=0,Zp2=0;   /*depth of pounded water on the surface*/
	double height_wind_sp;                        /*height of the Va measured for calculation of L */

	double Qhc_o[120],Qhc_u[120],Qhg[120]; /* The sensible heat flux from canopy and ground */
	double G[layer+2][120];                         /*   ' the heat flux into the canopy of over story --in W/m^2*/

	double Wcl_o[120],Wcs_o[120];   /* the masses od rain and snow on the canopy*/
	double Xcl_o[120],Xcs_o[120];  /* the fractoion of canopy coverd by liquid water and snow*/
	double Wcl_u[120],Wcs_u[120];   /* the masses of rain and snow on the canopy*/
	double Xcl_u[120],Xcs_u[120];   /* the fractoion of canopy coverd by liquid water and snow*/

	double r_rain_g[120]; /*the rainfall rate, on ground surface  m/s*/
	double rho_snow[120]; /*density of snow*/
	double alpha_v_sw[120],alpha_n_sw[120]; /*albeo of snow*/
	double Wg_snow[120];     /* the amount of snow on the ground*/
	double Xg_snow[120];     /* thw fraction of the ground surface overed by snow*/
	double Ac_snow_u[120];   /*the areas of canopy coverved by snow, o--overstory and  u-- understory*/
	double Ac_snow_o[120];

	double Ts0[120],Tsn0[120],Tsm0[120],Tsn1[120],Tsn2[120]; /*surface temperature*/
	double Tc_u[120]; /*the effective canopy temperature in K*/
	double Tm[layer+2][120]; /* Tb[layer+2][120],soil temeprature at the bottom and the middle of each layer*/

	double lambda[layer+2]; /*thermal conductivity of each soil layer;*/
	double Cs[layer+2][120];/*the soil voulumetric heat capacity of each soil layer, j+kkk/m^3/K*/

	double temp_air;    /* air temperature */
	double precip,rh_air,wind_sp;    /* precipitation in meters, relative humidity (%), windspeed in m/s */
	double temp_grd;    /* groud temperature */

	double Eil_o[120],EiS_o[120];/*the evaporation rate of intercepted water of overstory--in kg/m^2/s; l-- water; S-snow*/
	double Eil_u[120],EiS_u[120];/*the evaporation rate of intercepted water of overstory--in kg/m^2/s; l-- water; S-snowof intercepted water--in kg/m^2/s*/
	double Trans_o[120],Trans_u[120]; /*transpiration from overstory and understory canopies*/
	double Evap_soil[120];          /*evaporation from soil*/
	double Evap_SW[120];            /*evaporation from water pond*/
	double Evap_SS[120];            /*evaporation from snow pack*/
	
	double lambda_snow[120];   /*the effective thermal conductivity of snow --in m^2/s*/
	double e_a10;               /*the vapour partial pressure of water --in kPa(1mb=100Pa=0.1kpa)*/

	double Lv_liquid;                     /*the latent heat of vaporation from liquid at air temperature=Ta, in j+kkk/kg*/
	double Lv_solid=2.83*1000000;   /*the latent heat of vaporation from solid (snow/ice) at air temperature=Ta, in j+kkk/kg*/
	
	double  Ks;    /*KsCal,,KsMea[120]         instantaneous total short wave radiation (Global radiation)*/
	double  alpha_sat,alpha_dry;
	double  alpha_v_o,alpha_v_u;   /*visible albedo of overstory,  o--overstory, u--understory;*/
	double  alpha_n_o,alpha_n_u;   /*near-infrared albedoe of overstoryo,--overstory, u--understory;;*/
	double  alpha_g;                  /* the all-wave ground surface albedo*/
	double  alpha_v_g, alpha_n_g; /*the ground surface albedo for visble range and for near-infraded respectively*/
	double Cp_ca;                /*specific heat of moist air above the canopy  in j+kkk/k/G*/
	double ra_o,ra_u,ra_g;       /*the aerodynamic resistance of overstory, understory and groud surface   in s/m*/

	double  q_ca;             /*the actural canopy stomatal resistance  --in s/m*/
	double radiation_o, radiation_u, radiation_g;
	
	double ip=0;                            /*'the cumulative infiltration at the time of ponding   --in m/s*/
	double Inf=0;
	double zr=0.8;
	double Cpd=1004.65;
	double rho_w=1025.0;	// densoty of water

	double Cs_o_sunlit_old,Cs_o_shaded_old,Cs_u_sunlit_old,Cs_u_shaded_old; //CO2 concentration on the surfaces of leaves
	double Tc_o_sunlit_old,Tc_o_shaded_old,Tc_u_sunlit_old,Tc_u_shaded_old; //the effective canopy temperature in K

	double Gs_o_sunlit_new,Gs_o_shaded_new,Gs_u_sunlit_new,Gs_u_shaded_new; //stomatal conductance of the big leaf for water
	double Gs_o_sunlit_old,Gs_o_shaded_old,Gs_u_sunlit_old,Gs_u_shaded_old; //stomatal conductance of the big leaf for water

	double Ac_o_sunlit,Ac_o_shaded,Ac_u_sunlit,Ac_u_shaded;                 //net photosynthesis rate
	double Cs_o_sunlit_new,Cs_o_shaded_new,Cs_u_sunlit_new,Cs_u_shaded_new; //CO2 concentration on the surfaces of leaves
	double Ci_o_sunlit_new,Ci_o_shaded_new,Ci_u_sunlit_new,Ci_u_shaded_new; //intercellular CO2 concentration pn the leaf
	double Ci_o_sunlit_old,Ci_o_shaded_old,Ci_u_sunlit_old,Ci_u_shaded_old; //intercellular CO2 concentration pn the leaf
	double Cc_o_sunlit_new,Cc_o_shaded_new,Cc_u_sunlit_new,Cc_u_shaded_new; //CO2 concentration in the chloroplast

	double Tc_o_sunlit_new,Tc_o_shaded_new,Tc_u_sunlit_new,Tc_u_shaded_new; //the effective canopy temperature in K
	double f_soilwater;                                                          //an emperical parameter describing the relative availability of soil water for plants
	double Gw_o_sunlit,Gw_o_shaded,Gw_u_sunlit,Gw_u_shaded;                 //the total conductance for water from the intercellular space of the leaves to the reference height above the canopy
	double Gc_o_sunlit,Gc_o_shaded,Gc_u_sunlit,Gc_u_shaded;                 //the total conductance for CO2 from the intercellular space of the leaves to the reference height above the canopy
	double Gww_o_sunlit,Gww_o_shaded,Gww_u_sunlit,Gww_u_shaded;             //the total conductance for water from the surface of the leaves to the reference height above the canopy
	double Gh_o_sunlit,Gh_o_shaded,Gh_u_sunlit,Gh_u_shaded;                 //total conductance for heat transfer from the leaf surface to the reference height above the canopy

	double psychrometer;
	double R_o_sunlit,R_o_shaded,R_u_sunlit,R_u_shaded;                     //solar radiation absorbed by sunlit, shaded leaves

	double Tco, Tcu,slope;
	double H_o_sunlit,H_o_shaded;                     // sensible heat flux from leaves
	double LAI_o_sunlit,LAI_o_shaded,LAI_u_sunlit,LAI_u_shaded; //,stem_o,stem_u;
	double LAIo_sunlit,LAIo_shaded,LAIu_sunlit,LAIu_shaded;
	double radiation_o_sun, radiation_o_shaded;                             //net radiation of leaves
	double radiation_u_sun, radiation_u_shaded;
	double GPP_o_sunlit,GPP_o_shaded,GPP_u_sunlit,GPP_u_shaded;
        double SIF_o_sunlit,SIF_o_shaded,SIF_u_sunlit,SIF_u_shaded;             //SIF at canopy level   QBO
        double lSIF_o_sunlit,lSIF_o_shaded,lSIF_u_sunlit,lSIF_u_shaded;             //SIF at leaf level   QBO
        double stSIF_o_sunlit,stSIF_o_shaded,stSIF_u_sunlit,stSIF_u_shaded;             //SIF scat   QBO
        double AA,BB,CHIL,GMU,CVF;                                                       // QBO

	double VPS_air;
	double GH_o,G_o_a, G_o_b,G_u_a, G_u_b;
	double canopyh_o,canopyh_u;

	double VPD_air;     //Vapor pressure deficit term

	double mass_water_g;
    	double percentArea_snow_o, percentArea_snow_u;
	double Gheat_g;

	double b_h2o; // the intercept term in BWB model (mol H2O m-2 s-1)
	double m_h2o; // the slope in BWB model

	double leleaf_o_sunlit; // leaf latent heat flux (mol/m2/s)
	double leleaf_o_shaded; // leaf latent heat flux (mol/m2/s)
	double leleaf_u_sunlit; // leaf latent heat flux (mol/m2/s)
	double leleaf_u_shaded; // leaf latent heat flux (mol/m2/s)

	/* for the Vcmax - Nitrogen calculations  */
	double Kn = 0.3;  //0.713/2.4
	double G_theta = 0.5;
	double K,Vcmax0,Vcmax_sunlit,Vcmax_shaded,expr1,expr2,expr3;
	double slope_Vcmax_N, leaf_N,Jmax_sunlit,Jmax_shaded;

	psychrometer=0.066;

	alpha_sat = parameter[24];   // albedo of saturated /dry soil for module rainfall 1 
	alpha_dry = parameter[25];              /*the albedo of dry soil*/
	canopyh_o = parameter[29];	// to be used for module aerodynamic_conductance
	canopyh_u = parameter[30];
	height_wind_sp = parameter[31];	// the_height_to_measure_wind_speed, for module aerodynamic_conductance 
	m_h2o = parameter[33];    // to be used for module photosynthesis
	b_h2o = parameter[34];


	/* ---  Vcmax - Nitrogen calculations. by G.Mo  2011.04   ---  */

	if (CosZs>0)
	{
		K = G_theta*clumping/CosZs;
		Vcmax0 = parameter[36];
		expr1 = 1-exp(-K*lai);
		expr2 = 1-exp(-lai*(Kn+K));
		expr3 = 1-exp(-Kn*lai);

		if(expr1>0) Vcmax_sunlit = Vcmax0*parameter[47]*parameter[46]*K*expr2/(Kn+K)/expr1;
		else Vcmax_sunlit = Vcmax0;

		if (K>0 && lai>expr1/K) Vcmax_shaded = Vcmax0*parameter[47]*parameter[46]*(expr3/Kn-expr2/(Kn+K))/(lai-expr1/K);
		else Vcmax_shaded = Vcmax0;

	}
	
	
	/* ---  LAI calculation module by B. Chen   ---  */

	lai_o=lai;
	if (lai<0.1)  lai_o=0.1;
	landcover=(int) parameter[4];

	if (landcover==25 || landcover==40) lai_u=0.01;
	else	lai_u=1.18*exp(-0.99*lai_o);

	if (lai_u>lai_o) lai_u=0.01;

	stem_o=parameter[8]*0.2;
	stem_u=parameter[9]*0.2;

        /*    ------  G_fun ------    */
        CHIL=parameter[48];
        CVF=parameter[49];
        CVF=-1.0*log(lai)+CVF;
        AA = 0.5 - 0.633 * CHIL- 0.33 * CHIL* CHIL;
        BB = 0.877 * ( 1. - 2. * AA );
        GMU = ( AA + BB * CosZs ) / CosZs;

	lai2(stem_o,stem_u,landcover,CosZs,GMU,lai_o,clumping,lai_u,&LAIo_sunlit,&LAIo_shaded,&LAIu_sunlit,&LAIu_shaded,&LAI_o_sunlit,&LAI_o_shaded,&LAI_u_sunlit,&LAI_u_shaded); //Bin Chen



	/*  ---  initialization of this time step   --- */

	Ks = meteo->Srad;
	rh_air = meteo->rh; 
	wind_sp = meteo->wind;
	precip = meteo->rain/step;  /* precipitation in meters */
	temp_air=meteo->temp;
	
	if (Ks<=0)
	{
	 alpha_v_o = 0;	
	 alpha_n_o = 0;
	 alpha_v_u = 0;	
	 alpha_n_u = 0;
	}
	else
	{
	 alpha_v_o = parameter[22];
	 alpha_n_o = parameter[23];
	 alpha_v_u = parameter[22];
	 alpha_n_u = parameter[23];
	}


	Qhc_o[0]  = var_o[1] ;

	Ts0[0]  =var_o[3] ;        /* To:    The temperature of ground surface*/
	 if ((Ts0[0]-temp_air)>2.0) Ts0[0]=temp_air+2.0;
	 if ((Ts0[0]-temp_air)<-2.0) Ts0[0]=temp_air-2.0;
	Tsn0[0] =var_o[4] ;       /* To:    The temperature of ground surface*/
	 if ((Tsn0[0]-temp_air)>2.0) Tsn0[0]=temp_air+2.0;
	 if ((Tsn0[0]-temp_air)<-2.0) Tsn0[0]=temp_air-2.0;
	Tsm0[0] =var_o[5];
	 if ((Tsm0[0]-temp_air)>2.0) Tsm0[0]=temp_air+2.0;
	 if ((Tsm0[0]-temp_air)<-2.0) Tsm0[0]=temp_air-2.0;
	Tsn1[0] =var_o[6] ;       /* To:    The temperature of ground surface*/
	 if ((Tsn1[0]-temp_air)>2.0) Tsn1[0]=temp_air+2.0;
	 if ((Tsn1[0]-temp_air)<-2.0) Tsn1[0]=temp_air-2.0;
	Tsn2[0] =var_o[7] ;       /* To:    The temperature of ground surface*/
	 if ((Tsn2[0]-temp_air)>2.0) Tsn2[0]=temp_air+2.0;
	 if ((Tsn2[0]-temp_air)<-2.0) Tsn2[0]=temp_air-2.0;

	Wcl_o[0]=var_o[15] ;
	Wcs_o[0]=var_o[16] ;   /*the mass of intercepted liquid water and snow, overstory*/

	Wcl_u[0]=var_o[18] ;
	Wcs_u[0]=var_o[19] ;   /*the mass of intercepted liquid water and snow, overstory*/
	Wg_snow[0]=var_o[20]; /*thr fraction of ground surface covered by snow and snow mass*/

	soilp->Zsp = var_o[33];
	soilp->Zp = var_o[34];
	soilp->r_rain_g = var_o[35];
	
	Ac_snow_o[0] = var_o[36]; 
	Ac_snow_u[0] = var_o[37]; 
	rho_snow[0] = var_o[38]; 
	alpha_v_sw[0] = var_o[39]; 
	alpha_n_sw[0] = var_o[40]; 

	Zsp = soilp->Zsp;
	Zp = soilp->Zp;

    	if (Zp<0.001) Zp=0;

	for(i=9;i<=14;i++) soilp->temp_soil_p[i-9] = var_o[i];
	for(i=21;i<=26;i++) soilp->thetam_prev[i-21] = var_o[i];
	for(i=27;i<=32;i++) soilp->ice_ratio[i-27] = var_o[i];


	/* ---  vcmax jmax module  by L. He  --- */
	slope_Vcmax_N = parameter[47];
	leaf_N = parameter[46];

	Vcmax_Jmax(lai_o, clumping, Vcmax0,slope_Vcmax_N, leaf_N, CosZs, &Vcmax_sunlit, &Vcmax_shaded, &Jmax_sunlit, &Jmax_shaded);

	//temperatures of overstorey and understorey canopies
	Tc_o_sunlit_old=temp_air-0.5;
	Tc_o_shaded_old=temp_air-0.5;

	Tc_u_sunlit_old=temp_air-0.5;
	Tc_u_shaded_old=temp_air-0.5;


/*********************************************/
 	for(kkk=1;kkk<=kloop;kkk++)
 	{

	 /*  ---  Snow pack stage 1  by R.Luo  ---  */

    	snowpack_stage1 (temp_air, precip, Wcs_o[kkk-1],Wcs_u[kkk-1],Wg_snow[kkk-1],
                     &Wcs_o[kkk],&Wcs_u[kkk],&Wg_snow[kkk], lai_o,lai_u,clumping,
                     &Ac_snow_o[kkk], &Ac_snow_u[kkk], &Xcs_o[kkk], &Xcs_u[kkk], &Xg_snow[kkk] ,
                     &rho_snow[kkk], &Zsp, &alpha_v_sw[kkk], &alpha_n_sw[kkk]);

	/*  ---   rain fall stage 1   by R.Luo  ---   */
   
    	rainfall_stage1( temp_air,precip,Wcl_o[kkk-1],Wcl_u[kkk-1],
                    lai_o, lai_u, clumping, &Wcl_o[kkk], &Wcl_u[kkk], &Xcl_o[kkk],&Xcl_u[kkk],&r_rain_g[kkk]);


	 if(soilp->thetam_prev[1]<soilp->theta_vwp[1]*0.5) alpha_g = alpha_dry;
	 else alpha_g = (soilp->thetam_prev[1]-soilp->theta_vwp[1]*0.5)/(soilp->fei[1]-soilp->theta_vwp[1]*0.5) * (alpha_sat - alpha_dry) + alpha_dry;


	 alpha_v_g = 2.0 / 3.0 * alpha_g;     
	 alpha_n_g = 4.0 / 3.0 * alpha_g;


	 /*  --- soil_water_factor module by L. He ---  */
	 soil_water_factor_v2(soilp); 
	 f_soilwater = soilp->f_soilwater; 
	 
	 if (f_soilwater>1.0) f_soilwater=1.0; 	  // to be used for module photosynthesis

	 GH_o=Qhc_o[kkk-1];	// to be used as the init. for module aerodynamic_conductance

	 VPS_air=0.61078*exp(17.3*temp_air/(237.3+temp_air));  //to estimate saturated water vapor pressure in kpa
	 e_a10 =VPS_air*rh_air/100;		// to be used for module photosynthesis
	 VPD_air=VPS_air - e_a10;                                  //water vapor deficit at the reference height
	 
	 q_ca  =0.622 * e_a10 / (101.35- 0.378 * e_a10);                 //'in g/g  no dimention
	 Cp_ca =Cpd * (1 + 0.84 * q_ca);

	 slope       = 2503.0 / pow((temp_air+ 237.3),2) * exp(17.27 *temp_air/ (temp_air+ 237.3));

	 Gs_o_sunlit_old=1/200.0;Ci_o_sunlit_old=0.7*CO2_air;
	 Gs_o_shaded_old=1/200.0;Ci_o_shaded_old=0.7*CO2_air;
	 Gs_u_sunlit_old=1/300.0;Ci_u_sunlit_old=0.7*CO2_air;
	 Gs_u_shaded_old=1/300.0;Ci_u_shaded_old=0.7*CO2_air;

    	percentArea_snow_o=Ac_snow_o[kkk]/lai_o/2;
    	percentArea_snow_u=Ac_snow_u[kkk]/lai_u/2;
		 
	temp_grd = temp_air;   // ground temperature substituted by air temperature

	 num=0;
	 while(1)
	 {
		 num=num+1;

		 /* ---  aerodynamic_conductance module by G.Mo   ---  */

		 aerodynamic_conductance(canopyh_o,canopyh_u,height_wind_sp,clumping,temp_air,wind_sp,GH_o,lai_o+stem_o,lai_u+stem_u,&ra_o,&ra_u,&ra_g,
			 &G_o_a,&G_o_b,&G_u_a,&G_u_b);
		 

		 Gh_o_sunlit=1.0/(1.0/G_o_a+0.5/G_o_b);   //heat conductance of sunlit leaves of overstorey
		 Gh_o_shaded=1.0/(1.0/G_o_a+0.5/G_o_b);   //heat conductance of shaded leaves of overstorey

		 Gh_u_sunlit=1.0/(1.0/G_u_a+0.5/G_u_b);   //heat conductance of sunlit leaves of understorey
		 Gh_u_shaded=1.0/(1.0/G_u_a+0.5/G_u_b);   //heat conductance of shaded leaves of understorey

		 Gww_o_sunlit=1.0/(1.0/G_o_a+1.0/G_o_b+100); //conductance for intercepted water of sunlit leaves of overstorey
		 Gww_o_shaded=1.0/(1.0/G_o_a+1.0/G_o_b+100); //conductance for intercepted water of shaded leaves of overstorey
		 Gww_u_sunlit=1.0/(1.0/G_u_a+1.0/G_u_b+100); //conductance for intercepted water of sunlit leaves of understorey
		 Gww_u_shaded=1.0/(1.0/G_u_a+1.0/G_u_b+100); //conductance for intercepted water of shaded leaves of understorey

		 //temperatures of overstorey and understorey canopies
		 Tco=(Tc_o_sunlit_old*LAI_o_sunlit+Tc_o_shaded_old*LAI_o_shaded)/(LAI_o_sunlit+LAI_o_shaded);
		 Tcu=(Tc_u_sunlit_old*LAI_u_sunlit+Tc_u_shaded_old*LAI_u_shaded)/(LAI_u_sunlit+LAI_u_shaded);

 		/* ---  net Radiation at canopy and leaf level module by X.Luo   --- */


     		netRadiation (Ks,CosZs,Tco,Tcu,temp_grd,lai_o,lai_u,lai_o+stem_o,lai_u+stem_u,LAI_o_sunlit,LAI_o_shaded,LAI_u_sunlit,LAI_u_shaded,\
                   clumping,temp_air,rh_air,alpha_v_sw[kkk],alpha_n_sw[kkk],percentArea_snow_o,percentArea_snow_u, \
		   Xg_snow[kkk],alpha_v_o,alpha_n_o,alpha_v_u,alpha_n_u,alpha_v_g,alpha_n_g,\
                   &radiation_o,&radiation_u,&radiation_g,&radiation_o_sun, &radiation_o_shaded,&radiation_u_sun,&radiation_u_shaded,\
		   &R_o_sunlit,&R_o_shaded,&R_u_sunlit,&R_u_shaded);


 		/* ---  photosynthesis module by B. Chen   --- */
		 
		 Gw_o_sunlit = 1.0/(1.0/G_o_a+1.0/G_o_b+1.0/Gs_o_sunlit_old); //conductance of sunlit leaves of overstorey for water
	     	 Gw_o_shaded = 1.0/(1.0/G_o_a+1.0/G_o_b+1.0/Gs_o_shaded_old); //conductance of shaded leaves of overstorey for water
		 Gw_u_sunlit = 1.0/(1.0/G_u_a+1.0/G_u_b+1.0/Gs_u_sunlit_old); //conductance of sunlit leaves of understorey for water
		 Gw_u_shaded = 1.0/(1.0/G_u_a+1.0/G_u_b+1.0/Gs_u_shaded_old); //conductance of shaded leaves of understorey for water
		 
		 
		 leleaf_o_sunlit = Gw_o_sunlit*(VPD_air+slope*(Tc_o_sunlit_old-temp_air))*rho_a * Cp_ca/psychrometer; 
		 leleaf_o_shaded = Gw_o_shaded*(VPD_air+slope*(Tc_o_shaded_old-temp_air))*rho_a * Cp_ca/psychrometer;
		 leleaf_u_sunlit = Gw_u_sunlit*(VPD_air+slope*(Tc_u_sunlit_old-temp_air))*rho_a * Cp_ca/psychrometer;
		 leleaf_u_shaded = Gw_u_shaded*(VPD_air+slope*(Tc_u_shaded_old-temp_air))*rho_a * Cp_ca/psychrometer;

		 if (CosZs>0)
		 {
			photosynthesis(LAI_o_sunlit,Tc_o_sunlit_old,R_o_sunlit,e_a10,G_o_b,Vcmax_sunlit,f_soilwater,b_h2o,m_h2o,Ci_o_sunlit_old,temp_air,leleaf_o_sunlit,&Gs_o_sunlit_new,&Ac_o_sunlit,&Ci_o_sunlit_new,&lSIF_o_sunlit);
			photosynthesis(LAI_o_shaded,Tc_o_shaded_old,R_o_shaded,e_a10,G_o_b,Vcmax_shaded,f_soilwater,b_h2o,m_h2o,Ci_o_shaded_old,temp_air,leleaf_o_shaded,&Gs_o_shaded_new,&Ac_o_shaded,&Ci_o_shaded_new,&lSIF_o_shaded);
			photosynthesis(LAI_u_sunlit,Tc_u_sunlit_old,R_u_sunlit,e_a10,G_u_b,Vcmax_sunlit,f_soilwater,b_h2o,m_h2o,Ci_u_sunlit_old,temp_air,leleaf_u_sunlit,&Gs_u_sunlit_new,&Ac_u_sunlit,&Ci_u_sunlit_new,&lSIF_u_sunlit);
			photosynthesis(LAI_u_shaded,Tc_u_shaded_old,R_u_shaded,e_a10,G_u_b,Vcmax_shaded,f_soilwater,b_h2o,m_h2o,Ci_u_shaded_old,temp_air,leleaf_u_shaded,&Gs_u_shaded_new,&Ac_u_shaded,&Ci_u_shaded_new,&lSIF_u_shaded);
		 }
		 else
		 {
			 Gs_o_sunlit_new=0.0001;
			 Ac_o_sunlit=0.0;
			 Ci_o_sunlit_new=CO2_air*0.7;
			 Cs_o_sunlit_new=CO2_air;
			 Cc_o_sunlit_new=CO2_air*0.7*0.8;
                         lSIF_o_sunlit=0.0;

			 Gs_o_shaded_new=0.0001;
			 Ac_o_shaded=0.0;
			 Ci_o_shaded_new=CO2_air*0.7;
			 Cs_o_shaded_new=CO2_air;
			 Cc_o_shaded_new=CO2_air*0.7*0.8;
                         lSIF_o_shaded=0.0;

			 Gs_u_sunlit_new=0.0001;
			 Ac_u_sunlit=0.0;
			 Ci_u_sunlit_new=CO2_air*0.7;
			 Cs_u_sunlit_new=CO2_air;
			 Cc_u_sunlit_new=CO2_air*0.7*0.8;
                         lSIF_u_sunlit=0.0;

			 Gs_u_shaded_new=0.0001;
			 Ac_u_shaded=0.0;
			 Ci_u_shaded_new=CO2_air*0.7;
			 Cs_u_shaded_new=CO2_air;
			 Cc_u_shaded_new=CO2_air*0.7*0.8;
                         lSIF_u_shaded=0.0;

		 }

		 Ci_o_sunlit_old=Ci_o_sunlit_new;
		 Cs_o_sunlit_old=Cs_o_sunlit_new;
		 Gs_o_sunlit_old=Gs_o_sunlit_new; // m/s
		 Gw_o_sunlit=1.0/(1.0/G_o_a+1.0/G_o_b+1.0/Gs_o_sunlit_new); //conductance of sunlit leaves of overstorey for wwater
		 Gc_o_sunlit=1.0/(1.0/G_o_a+1.4/G_o_b+1.6/Gs_o_sunlit_new); //conductance of sunlit leaves of overstorey for CO2
			 
	     	 Ci_o_shaded_old=Ci_o_shaded_new;
		 Cs_o_shaded_old=Cs_o_shaded_new;
		 Gs_o_shaded_old=Gs_o_shaded_new;   // m/s
		 Gw_o_shaded=1.0/(1.0/G_o_a+1.0/G_o_b+1.0/Gs_o_shaded_new); //conductance of sunlit leaves of overstorey for wwater
		 Gc_o_shaded=1.0/(1.0/G_o_a+1.4/G_o_b+1.6/Gs_o_shaded_new); //conductance of sunlit leaves of overstorey for CO2


		 Ci_u_sunlit_old=Ci_o_sunlit_new;
		 Cs_u_sunlit_old=Cs_u_sunlit_new;
		 Gs_u_sunlit_old=Gs_o_sunlit_new;   // m/s
		 Gw_u_sunlit=1.0/(1.0/G_u_a+1.0/G_u_b+1.0/Gs_u_sunlit_new); //conductance of sunlit leaves of overstorey for wwater
		 Gc_u_sunlit=1.0/(1.0/G_u_a+1.4/G_u_b+1.6/Gs_u_sunlit_new); //conductance of sunlit leaves of overstorey for CO2


		 Ci_u_shaded_old=Ci_u_shaded_new;
		 Cs_u_shaded_old=Cs_u_shaded_new;
		 Gs_u_shaded_old=Gs_u_shaded_new;   // m/s
		 Gw_u_shaded=1.0/(1.0/G_u_a+1.0/G_u_b+1.0/Gs_u_shaded_new); //conductance of sunlit leaves of overstorey for wwater
		 Gc_u_shaded=1.0/(1.0/G_u_a+1.4/G_u_b+1.6/Gs_u_shaded_new); //conductance of sunlit leaves of overstorey for CO2



		/* --- leaf tempratures module by L. He  --- */

		Leaf_Temperatures(temp_air, slope, psychrometer, VPD_air, Cp_ca, Gw_o_sunlit, Gw_o_shaded, Gw_u_sunlit, Gw_u_shaded, \
		  Gww_o_sunlit, Gww_o_shaded, Gww_u_sunlit, Gww_u_shaded, Gh_o_sunlit, Gh_o_shaded, Gh_u_sunlit, Gh_u_shaded, \
		  Xcs_o[kkk], Xcl_o[kkk], Xcs_u[kkk], Xcl_u[kkk], radiation_o_sun, radiation_o_shaded, radiation_u_sun, radiation_u_shaded, \
		  &Tc_o_sunlit_new,&Tc_o_shaded_new,&Tc_u_sunlit_new,&Tc_u_shaded_new);


		 H_o_sunlit =(Tc_o_sunlit_new-temp_air)*rho_a * Cp_ca*Gh_o_sunlit; 
		 H_o_shaded =(Tc_o_shaded_new-temp_air)*rho_a * Cp_ca*Gh_o_shaded;
		 GH_o=H_o_sunlit*LAI_o_sunlit+H_o_shaded*LAI_o_shaded; 	// for next num aerodynamic conductance calculation

		 if(fabs(Tc_o_sunlit_new-Tc_o_sunlit_old)<0.02 && fabs(Tc_o_shaded_new-Tc_o_shaded_old)<0.02 &&
			 fabs(Tc_u_sunlit_new-Tc_u_sunlit_old)<0.02 && fabs(Tc_u_shaded_new-Tc_u_shaded_old)<0.02 )
 
			break;
		 else
			 if(num>22)  //if the iteration does not converge.
			 {

				 Tc_o_sunlit_old=temp_air;
				 Tc_o_shaded_old=temp_air;
				 Tc_u_sunlit_old=temp_air;
				 Tc_u_shaded_old=temp_air;
				 break;
			 }
			 else
			 {
				 Tc_o_sunlit_old=Tc_o_sunlit_new;
				 Tc_o_shaded_old=Tc_o_shaded_new;
				 Tc_u_sunlit_old=Tc_u_sunlit_new;
				 Tc_u_shaded_old=Tc_u_shaded_new;
			 }

	 }	// end of while

	 GPP_o_sunlit=Ac_o_sunlit*LAIo_sunlit;
	 GPP_o_shaded=Ac_o_shaded*LAIo_shaded;
	 GPP_u_sunlit=Ac_u_sunlit*LAIu_sunlit;
	 GPP_u_shaded=Ac_u_shaded*LAIu_shaded;


         stSIF_o_sunlit=lSIF_o_sunlit*0.5*clumping*(1.1-0.1*LAIo_sunlit)*exp(-CosZs)*LAIo_sunlit;
         stSIF_o_shaded=lSIF_o_shaded*0.5*clumping*(1.1-0.1*LAIo_shaded)*exp(-CosZs)*LAIo_shaded;
         stSIF_u_sunlit=lSIF_u_sunlit*0.5*clumping*(1.1-0.1*LAIu_sunlit)*exp(-CosZs)*LAIu_sunlit;
         stSIF_u_shaded=lSIF_u_shaded*0.5*clumping*(1.1-0.1*LAIu_shaded)*exp(-CosZs)*LAIu_shaded;

         SIF_o_sunlit=lSIF_o_sunlit*exp(-GMU*clumping*LAIo_sunlit/CosZs)*LAIo_sunlit+stSIF_o_sunlit;
         SIF_o_shaded=lSIF_o_shaded*exp(-GMU*clumping*LAIo_shaded/CosZs)*LAIo_shaded+stSIF_o_shaded;
         SIF_u_sunlit=lSIF_u_sunlit*exp(-GMU*clumping*LAIu_sunlit/CosZs)*LAIu_sunlit+stSIF_u_sunlit;
         SIF_u_shaded=lSIF_u_shaded*exp(-GMU*clumping*LAIu_shaded/CosZs)*LAIu_shaded+stSIF_u_shaded;


	 /*  ---   Transpiration module by X. Luo   ---  */

	 transpiration(Tc_o_sunlit_new, Tc_o_shaded_new, Tc_u_sunlit_new, Tc_u_shaded_new,temp_air, rh_air,
                    Gw_o_sunlit, Gw_o_shaded, Gw_u_sunlit, Gw_u_shaded,LAI_o_sunlit, LAI_o_shaded, LAI_u_sunlit, LAI_u_shaded,
                    &Trans_o[kkk], &Trans_u[kkk]);

	 /*  ---   Evaporation and sublimation from canopy by X. Luo  ---   */

	 evaporation_canopy (Tc_o_sunlit_new, Tc_o_shaded_new, Tc_u_sunlit_new, Tc_u_shaded_new,temp_air, rh_air,
                    Gww_o_sunlit, Gww_o_shaded, Gww_u_sunlit, Gww_u_shaded,
                    LAI_o_sunlit, LAI_o_shaded, LAI_u_sunlit, LAI_u_shaded,
                    Xcl_o[kkk], Xcl_u[kkk], Xcs_o[kkk], Xcs_u[kkk],
                    &Eil_o[kkk], &Eil_u[kkk], &EiS_o[kkk], &EiS_u[kkk]
                    );

	 /*  ---   Rainfall stage 2 by X. Luo  ---  */

	 rainfall_stage2(Eil_o[kkk], Eil_u[kkk], &Wcl_o[kkk], &Wcl_u[kkk]);


	 /*  ---  Snow pack stage 2  by X. Luo  ---  */

	 snowpack_stage2(EiS_o[kkk], EiS_u[kkk], &Wcs_o[kkk], &Wcs_u[kkk]);

	 
    	 /*  ---  Evaporation from soil module by X. Luo  ---  */
    	 Gheat_g=1/ra_g;
    	 mass_water_g=rho_w*Zp;

    	 evaporation_soil(temp_grd, Ts0[kkk-1], rh_air, radiation_g, Gheat_g, &Xg_snow[kkk],
                       &Zp, &Zsp, &mass_water_g, &Wg_snow[kkk], rho_snow[kkk], soilp->thetam_prev[0], soilp->fei[0],
                       &Evap_soil[kkk], &Evap_SW[kkk],&Evap_SS[kkk]);


/*  to be checked later:  why set these 4 to 0?  */
	 Eil_o[kkk]=0;EiS_o[kkk]=0;
	 Eil_u[kkk]=0;EiS_u[kkk]=0;

	 /*  --- soil Thermal Conductivity module by L. He ---  */

	 UpdateSoilThermalConductivity(soilp); 

	 Update_Cs(soilp); 


	 /*  ---   Surface temperature by X. Luo  ---  */

	 Cs[0][kkk]=soilp->Cs[0];  //added
	 Cs[1][kkk]=soilp->Cs[0];
	 Tc_u[kkk]=Tcu;		//added
	 lambda[1]=soilp->lambda[0];
	 d_soil[1]=soilp->d_soil[0];
	 Tm[1][kkk-1]=soilp->temp_soil_p[1]; // first place is temp_soil_p[0]?
	 Tm[0][kkk-1]=soilp->temp_soil_p[0];
	 G[1][kkk]=soilp->G[0];

    	 surface_temperature (temp_air,rh_air,Zsp,Zp,
                          Cs[1][kkk], Cs[0][kkk], Gheat_g, d_soil[1], rho_snow[kkk],Tc_u[kkk],
                          radiation_g,Evap_soil[kkk],Evap_SW[kkk],Evap_SS[kkk],
                          lambda[1],Xg_snow[kkk],G[1][kkk],
                          Ts0[kkk-1],Tm[1][kkk-1],Tm[0][kkk-1], Tsn0[kkk-1],
                          Tsm0[kkk-1],Tsn1[kkk-1], Tsn2[kkk-1],
                          &Ts0[kkk], &Tm[0][kkk], &Tsn0[kkk],
                          &Tsm0[kkk], &Tsn1[kkk], &Tsn2[kkk],
                          &G[0][kkk]);

    	 soilp->temp_soil_c[0]=Tm[0][kkk];

	

	 /*  ---   Snow Pack Stage 3 module by X. Luo   ---  */

    	 snowpack_stage3(temp_air, Tsn0[kkk], Tsn0[kkk-1],rho_snow[kkk], &Zsp, &Zp, &Wg_snow[kkk]);


	 /*  ---  Sensible heat flux module by X. Luo    ---  */

    	 sensible_heat(Tc_o_sunlit_new, Tc_o_shaded_new, Tc_u_sunlit_new, Tc_u_shaded_new, Ts0[kkk], temp_air, rh_air,
                    Gh_o_sunlit, Gh_o_shaded, Gh_u_sunlit, Gh_u_shaded, Gheat_g,
                    LAI_o_sunlit, LAI_o_shaded, LAI_u_sunlit, LAI_u_shaded,
                    &Qhc_o[kkk], &Qhc_u[kkk], &Qhg[kkk]);

	 
	 /*  --- soil water module by L. He ---  */

	 soilp->Zsp = Zsp;
	 soilp->G[0] = G[0][kkk];

	 UpdateHeatFlux(soilp, Xg_snow[kkk], lambda_snow[kkk], Tsn0[kkk], temp_air, kstep); 

	 Soil_Water_Uptake(soilp, Trans_o[kkk], Trans_u[kkk], Evap_soil[kkk]); 

	 soilp->r_rain_g = r_rain_g[kkk];
	 soilp->Zp = Zp;

	 UpdateSoilMoisture(soilp, kstep);

	 Zp = soilp->Zp;

  	}   /* The end of kkk loop  */


  	kkk=kloop ;

  	if (Tsn1[kkk]>40) Tsn2[kkk]=40;
  	if (Tsn1[kkk]<-40) Tsn2[kkk]=-40;

  	if (Tsn2[kkk]>40) Tsn2[kkk]=40;
  	if (Tsn2[kkk]<-40) Tsn2[kkk]=-40;

  	var_n[1]=Qhc_o[kkk];	// sensible heta flux to be used by next step
 
  	var_n[3]=Ts0[kkk];        /* To:    The temperature of ground surface*/
  	var_n[4]=Tsn0[kkk];       /* To:    The temperature of ground surface*/
  	var_n[5]=Tsm0[kkk];
  	var_n[6]=Tsn1[kkk];       /* To:    The temperature of ground surface*/
  	var_n[7]=Tsn2[kkk];       /* To:    The temperature of ground surface*/

 	for(i=9;i<=14;i++)  var_n[i] = 	soilp->temp_soil_c[i-9];
	for(i=21;i<=26;i++) var_n[i] = soilp->thetam[i-21];
	for(i=27;i<=32;i++) var_n[i] = soilp->ice_ratio[i-27];


  	var_n[15]=Wcl_o[kkk];
  	var_n[16]=Wcs_o[kkk];   /*the mass of intercepted liquid water and snow, overstory*/

  	var_n[18]=Wcl_u[kkk];
  	var_n[19]=Wcs_u[kkk];   /*the mass of intercepted liquid water and snow, overstory*/
  	var_n[20]=Wg_snow[kkk]; /*thr fraction of ground surface covered by snow and snow mass*/

	var_n[33] = soilp->Zsp;
	var_n[34] = soilp->Zp;
	var_n[35] = soilp->r_rain_g;

	var_n[36] = Ac_snow_o[kkk]; 
	var_n[37] = Ac_snow_u[kkk]; 
	var_n[38] = rho_snow[kkk]; 
	var_n[39] = alpha_v_sw[kkk]; 
	var_n[40] = alpha_n_sw[kkk]; 


	Lv_liquid = (2.501 - 0.00237 * temp_air) * 1000000;  /*The laten heat of water vaporization in j/kg */

  	mid_res->Net_Rad = radiation_o+radiation_u+radiation_g;

  	mid_res->LH = Lv_liquid *(Trans_o[kkk]+Eil_o[kkk]+Trans_u[kkk]+Eil_u[kkk]+Evap_soil[kkk]+Evap_SW[kkk])+Lv_solid*(EiS_o[kkk]+EiS_u[kkk]+Evap_SS[kkk]);
     	//  LH:  total latent heat flux
  	
	mid_res->SH = Qhc_o[kkk]+Qhc_u[kkk]+Qhg[kkk]; // SH: total sensible heta flux 

  	mid_res->Trans = (Trans_o[kkk]+Trans_u[kkk])*step;	// total transpiration  mm/step
  	mid_res->Evap = (Eil_o[kkk]+Eil_u[kkk]+Evap_soil[kkk]+Evap_SW[kkk]+EiS_o[kkk]+EiS_u[kkk]+Evap_SS[kkk])*step;   // total evapiration  mm/step 

	mid_res->gpp_o_sunlit = GPP_o_sunlit;   // umol C/m2/s
	mid_res->gpp_u_sunlit = GPP_u_sunlit; 
	mid_res->gpp_o_shaded = GPP_o_shaded;	
	mid_res->gpp_u_shaded = GPP_u_shaded; 

  	// total GPP gC/m2/step
  	mid_res->GPP = (GPP_o_sunlit+GPP_o_shaded+GPP_u_sunlit+GPP_u_shaded)*12*step*0.000001;
        // total SIF
        mid_res->SIF = (SIF_o_sunlit+SIF_o_shaded+SIF_u_sunlit+SIF_u_shaded)/CVF;

//  if (rstep==1)  printf("t=%f s=%f\n",var_n[10],var_n[22]);


return;

}
