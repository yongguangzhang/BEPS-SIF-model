// Header file for struct Soil
// Created on Dec. 05, 2012 by Liming He
// Last revision on May 15, 2015

#ifndef SOIL_H
#define SOIL_H

#define FW_VERSION 1 // 0 for soil water uptake using R, and 1 for soil water uptake using R*fpsisr

#define max(a,b) 	((a)>(b))?(a):(b) 
#define min(a,b)	((a)<(b))?(a):(b) 

// Note: change the value of MAX_LAYERS to a small one for global application
// e. g. max layers = 6.
#define MAX_LAYERS 10 // LHE. Jan 28, 2013.
#define DEPTH_F 6 

struct Soil{
	//////////////////////////////////////////////////////////////////////////////////////////////////////
	// Properties belong to the whole soil profile
	int flag;				// reserved for EnKF usage. 
	int n_layer;				/* the number of layers used in the model. Make sure n_layer <= MAX_LAYERS */
	int step_peroid;

	/* condition on the top boundary */
	double Zp;				 /* depth of ponded water on the ground surface*/
	double Zsp;              /*Snow depth;*/
	double r_rain_g;		 /*the rainfall rate, un--on understory g--on ground surface  m/s */
	double soil_r;			// soil surface resistance for water, disuss with Remi - an interface here. 
	double r_drainage;

	/* some variable used for soil */
	//double t1;
	//double t2;
	double r_root_decay;			/* decay_rate_of_root_distribution */
	double psi_min; // for fw
	double alpha; // for fw
	double f_soilwater;

	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Properties belong to each soil horizon
	double d_soil[MAX_LAYERS];
	double f_root[MAX_LAYERS]; /* root weight*/
	double dt[MAX_LAYERS]; /* the weight calcualted from soil_water_factor **re-calculate in the model*/
	// from read-param
	double thermal_cond[MAX_LAYERS]; /* thermal conductivity. Unit:   */
	double theta_vfc[MAX_LAYERS]; /* field capacity */ // not used in this model. LHE. Feb. 01, 2013
	double theta_vwp[MAX_LAYERS]; /* wilt point*/
	double fei[MAX_LAYERS];  /* porosity */
	double Ksat[MAX_LAYERS]; /* saturated hydraulic conductivity */
	double psi_sat[MAX_LAYERS]; /* water potential at sat */
	double b[MAX_LAYERS]; /* Cambell parameter b */
	double density_soil[MAX_LAYERS]; // soil bulk density of layer. LHE. Feb. 12, 2013.
	double f_org[MAX_LAYERS]; // volume fraction of organic matter in layer (%). 

	// need to save
	double ice_ratio[MAX_LAYERS];/* The ratio of ice of soil layer **save */
	double thetam[MAX_LAYERS], thetam_prev[MAX_LAYERS]; /* soil water content in this layer */
	double temp_soil_p[MAX_LAYERS]; /* soil temperature in this layer */ // don't change because it is used in soil_water_factor_v2, and UpdateSoil_Moisture.
	double temp_soil_c[MAX_LAYERS]; /* soil temperature in this layer */ // don't change because it is used in soil_water_factor_v2, and UpdateSoil_Moisture.

	// derived variables below:
	double f_ice[MAX_LAYERS]; // derived var.
	double psim[MAX_LAYERS]; /*soil water suction in this layer. Note: this variable can be derived from other parameters. LHE. */
	double thetab[MAX_LAYERS]; /* soil water content at the bottom of each layer*/
	double psib[MAX_LAYERS]; /*soil water suction at the bottom this layer.  */
	double r_waterflow[MAX_LAYERS]; /*'the liquid water flow rates at the soil layer interfaces  'eg. 0,1,2..., represents the surface, the bottom of layer1, the bottom of layer2,...*/

	double km[MAX_LAYERS], Kb[MAX_LAYERS];   /*'the hydraulic conductivity of certain soil layer*/
	double KK[MAX_LAYERS];           /* The average  conductivity of two soil layers.*/ 

	double Cs[MAX_LAYERS];
	double lambda[MAX_LAYERS];  /*thermal conductivity of each soil layer;*/ /* ={0} by LHE */ // not used in gpp-only version. derived var.
	double Ett[MAX_LAYERS]; /*ET in each layer. derived var. */
	/* define a lambda_top for ice? */
	double G[MAX_LAYERS]; /* energy fluxes */
};

void SoilRootFraction(struct Soil soil[]);
void Init_Soil_Parameters(int landcover, int stxt, double r_root_decay, struct Soil p[]);
void Init_Soil_Status(struct Soil p[], double Tsoil, double Tair, double Ms, double snowdepth);
void soil_water_factor_v2(struct Soil p[]);
void Soil_Water_Uptake(struct Soil p[], double Trans_o, double Trans_u, double Evap_soil);
void UpdateSoilLambda(struct Soil soil[]);
void init_soil_parameter(unsigned char T_USDA, unsigned char S_USDA, unsigned char Ref_Depth, double T_Density, double S_Density, double T_OC, double S_OC, struct Soil soil[]);
void Update_Cs(struct Soil p[]);
void Update_ice_ratio(struct Soil p[]);
void UpdateSoilThermalConductivity(struct Soil p[]);
void UpdateHeatFlux(struct Soil p[], double Xg_snow, double lambda_snow, double Tsn0, double Tair_annual_mean, double peroid_in_seconds);
void UpdateSoilMoisture(struct Soil p[], double peroid_in_seconds);
#endif

