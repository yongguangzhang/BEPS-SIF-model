#include "soil.h"
#include <math.h>

// Liming He
// June 2, 2015

// function: 
//			(1) set the depth for each layer
//			(2) set the parameters for each layer

void Init_Soil_Parameters(int landcover, int stxt, double r_root_decay, struct Soil p[])
{
	p->n_layer = 5;

	if (landcover == 6 || landcover == 9)
	{
		p->psi_min = 10.0;
		p->alpha = 1.5;
	}
	else
	{
		p->psi_min = 33.0;
		p->alpha = 0.4;
	}

	p->d_soil[0] = 0.05;		/* depth_layer0 */
	p->d_soil[1] = 0.10;		/* depth_layer1 */
	p->d_soil[2] = 0.20;		/* depth_layer2 */
	p->d_soil[3] = 0.40;		/* depth_layer3 */
	p->d_soil[4] = 1.25;		/* depth_layer4 */

	p->r_root_decay = r_root_decay;
	SoilRootFraction(p);
	
	p->density_soil[0] = 1300.0;   // from flux tower.
	p->density_soil[1] = 1500.0;
	p->density_soil[2] = 1517.0;
	p->density_soil[3] = 1517.0;
	p->density_soil[4] = 1517.0;
	
	p->f_org[0] = 5;
	p->f_org[1] = 2;
	p->f_org[2] = 1;
	p->f_org[3] = 1;
	p->f_org[4] = 0.3;
	


	switch (stxt)
	{
	case 1:  //sand
		p->b[0] = 1.7;          p->b[1] = 1.9;           p->b[2] = 2.1;             p->b[3] = 2.3;          p->b[4] = 2.5;
		p->Ksat[0] = 0.000058;      p->Ksat[1] = 0.000052;       p->Ksat[2] = 0.000046;         p->Ksat[3] = 0.000035;      p->Ksat[4] = 0.000010;       /* saturated hydraulic conductivity */
		p->fei[0] = 0.437;      p->fei[1] = 0.437;       p->fei[2] = 0.437;         p->fei[3] = 0.437;      p->fei[4] = 0.437;       /* porosity */
		p->theta_vfc[0] = 0.09; p->theta_vfc[1] = 0.09;  p->theta_vfc[2] = 0.09;    p->theta_vfc[3] = 0.09; p->theta_vfc[4] = 0.09;   /* field capacity */
		p->theta_vwp[0] = 0.03; p->theta_vwp[1] = 0.03;  p->theta_vwp[2] = 0.03;    p->theta_vwp[3] = 0.03; p->theta_vwp[4] = 0.03;  /* wilt point*/
		p->thermal_cond[0] = 8.6;  p->thermal_cond[1] = 8.6;   p->thermal_cond[2] = 8.6;     p->thermal_cond[3] = 8.6;  p->thermal_cond[4] = 8.6;  /* thermal conductivity */
		p->psi_sat[0] = 0.07;   p->psi_sat[1] = 0.08;    p->psi_sat[2] = 0.09;      p->psi_sat[3] = 0.10;   p->psi_sat[4] = 0.12;  /* water potential at sat */
		break;

	case 2: //loamy sand
		p->b[0] = 2.1;          p->b[1] = 2.3;           p->b[2] = 2.5;             p->b[3] = 2.7;          p->b[4] = 2.9;
		p->Ksat[0] = 0.000017;      p->Ksat[1] = 0.000015;       p->Ksat[2] = 0.000014;         p->Ksat[3] = 0.000010;      p->Ksat[4] = 0.000003;
		p->fei[0] = 0.437;      p->fei[1] = 0.437;       p->fei[2] = 0.437;         p->fei[3] = 0.437;      p->fei[4] = 0.437;
		p->theta_vfc[0] = 0.21; p->theta_vfc[1] = 0.21;  p->theta_vfc[2] = 0.21;    p->theta_vfc[3] = 0.21; p->theta_vfc[4] = 0.21;
		p->theta_vwp[0] = 0.06; p->theta_vwp[1] = 0.06;  p->theta_vwp[2] = 0.06;    p->theta_vwp[3] = 0.06; p->theta_vwp[4] = 0.06;
		p->thermal_cond[0] = 8.3;  p->thermal_cond[1] = 8.3;   p->thermal_cond[2] = 8.3;     p->thermal_cond[3] = 8.3;  p->thermal_cond[4] = 8.3;
		p->psi_sat[0] = 0.09;   p->psi_sat[1] = 0.10;    p->psi_sat[2] = 0.11;      p->psi_sat[3] = 0.12;   p->psi_sat[4] = 0.14;
		break;

	case 3: // sandy loam
		p->b[0] = 3.1;          p->b[1] = 3.3;           p->b[2] = 3.5;             p->b[3] = 3.7;          p->b[4] = 3.9;
		p->Ksat[0] = 0.0000072;     p->Ksat[1] = 0.00000648;     p->Ksat[2] = 0.00000576;       p->Ksat[3] = 0.00000432;    p->Ksat[4] = 0.00000144;
		p->fei[0] = 0.453;      p->fei[1] = 0.453;       p->fei[2] = 0.453;         p->fei[3] = 0.453;      p->fei[4] = 0.453;
		p->theta_vfc[0] = 0.21; p->theta_vfc[1] = 0.21;  p->theta_vfc[2] = 0.21;    p->theta_vfc[3] = 0.21; p->theta_vfc[4] = 0.21;
		p->theta_vwp[0] = 0.10; p->theta_vwp[1] = 0.10;  p->theta_vwp[2] = 0.10;    p->theta_vwp[3] = 0.10; p->theta_vwp[4] = 0.10;
		p->thermal_cond[0] = 8.0;  p->thermal_cond[1] = 8.0;   p->thermal_cond[2] = 8.0;     p->thermal_cond[3] = 8.0;  p->thermal_cond[4] = 8.0;
		p->psi_sat[0] = 0.15;   p->psi_sat[1] = 0.16;    p->psi_sat[2] = 0.17;      p->psi_sat[3] = 0.18;   p->psi_sat[4] = 0.20;
		break;

	case 4: // loam
		p->b[0] = 4.5;          p->b[1] = 4.7;           p->b[2] = 4.9;             p->b[3] = 5.1;          p->b[4] = 5.3;
		p->Ksat[0] = 0.0000037;     p->Ksat[1] = 0.0000033;      p->Ksat[2] = 0.00000296;       p->Ksat[3] = 0.00000222;    p->Ksat[4] = 0.00000074;
		p->fei[0] = 0.463;    p->fei[1] = 0.463;       p->fei[2] = 0.463;         p->fei[3] = 0.463;      p->fei[4] = 0.463;
		p->theta_vfc[0] = 0.27; p->theta_vfc[1] = 0.27;  p->theta_vfc[2] = 0.27;    p->theta_vfc[3] = 0.27; p->theta_vfc[4] = 0.27;
		p->theta_vwp[0] = 0.12; p->theta_vwp[1] = 0.12;  p->theta_vwp[2] = 0.12;    p->theta_vwp[3] = 0.12; p->theta_vwp[4] = 0.12;
		p->thermal_cond[0] = 7.0;  p->thermal_cond[1] = 7.0;   p->thermal_cond[2] = 7.0;     p->thermal_cond[3] = 7.0;  p->thermal_cond[4] = 7.0;
		p->psi_sat[0] = 0.11;   p->psi_sat[1] = 0.12;    p->psi_sat[2] = 0.13;      p->psi_sat[3] = 0.14;   p->psi_sat[4] = 0.16;
		break;

	case 5: // silty loam
		p->b[0] = 4.7;          p->b[1] = 4.9;           p->b[2] = 5.1;             p->b[3] = 5.3;          p->b[4] = 5.5;
		p->Ksat[0] = 0.0000019;     p->Ksat[1] = 0.0000017;      p->Ksat[2] = 0.00000152;       p->Ksat[3] = 0.00000114;    p->Ksat[4] = 0.00000038;
		p->fei[0] = 0.501;      p->fei[1] = 0.501;       p->fei[2] = 0.501;         p->fei[3] = 0.501;      p->fei[4] = 0.501;
		p->theta_vfc[0] = 0.33; p->theta_vfc[1] = 0.33;  p->theta_vfc[2] = 0.33;    p->theta_vfc[3] = 0.33; p->theta_vfc[4] = 0.33;
		p->theta_vwp[0] = 0.13; p->theta_vwp[1] = 0.13;  p->theta_vwp[2] = 0.13;    p->theta_vwp[3] = 0.13; p->theta_vwp[4] = 0.13;
		p->thermal_cond[0] = 6.3;  p->thermal_cond[1] = 6.3;   p->thermal_cond[2] = 6.3;     p->thermal_cond[3] = 6.3;  p->thermal_cond[4] = 6.3;
		p->psi_sat[0] = 0.21;   p->psi_sat[1] = 0.22;    p->psi_sat[2] = 0.23;      p->psi_sat[3] = 0.24;   p->psi_sat[4] = 0.26;
		break;

	case 6: // sandy caly loam
		p->b[0] = 4.0;          p->b[1] = 4.2;           p->b[2] = 4.4;             p->b[3] = 4.6;          p->b[4] = 4.8;
		p->Ksat[0] = 0.0000012;     p->Ksat[1] = 0.00000108;     p->Ksat[2] = 0.0000096;        p->Ksat[3] = 0.0000072;     p->Ksat[4] = 0.0000024;
		p->fei[0] = 0.398;      p->fei[1] = 0.398;       p->fei[2] = 0.398;         p->fei[3] = 0.398;      p->fei[4] = 0.398;
		p->theta_vfc[0] = 0.26; p->theta_vfc[1] = 0.26;  p->theta_vfc[2] = 0.26;    p->theta_vfc[3] = 0.26; p->theta_vfc[4] = 0.26;
		p->theta_vwp[0] = 0.15; p->theta_vwp[1] = 0.15;  p->theta_vwp[2] = 0.15;    p->theta_vwp[3] = 0.15; p->theta_vwp[4] = 0.15;
		p->thermal_cond[0] = 7.0;  p->thermal_cond[1] = 7.0;   p->thermal_cond[2] = 7.0;     p->thermal_cond[3] = 7.0;  p->thermal_cond[4] = 7.0;
		p->psi_sat[0] = 0.28;   p->psi_sat[1] = 0.29;    p->psi_sat[2] = 0.30;      p->psi_sat[3] = 0.31;   p->psi_sat[4] = 0.33;
		break;


	case 7: // clay loam
		p->b[0] = 5.2;          p->b[1] = 5.4;           p->b[2] = 5.6;             p->b[3] = 5.8;          p->b[4] = 6.0;
		p->Ksat[0] = 0.00000064;    p->Ksat[1] = 0.00000058;     p->Ksat[2] = 0.00000051;       p->Ksat[3] = 0.00000038;    p->Ksat[4] = 0.00000013;
		p->fei[0] = 0.464;      p->fei[1] = 0.464;       p->fei[2] = 0.464;         p->fei[3] = 0.464;      p->fei[4] = 0.464;
		p->theta_vfc[0] = 0.32; p->theta_vfc[1] = 0.32;  p->theta_vfc[2] = 0.32;    p->theta_vfc[3] = 0.32; p->theta_vfc[4] = 0.32;
		p->theta_vwp[0] = 0.20; p->theta_vwp[1] = 0.20;  p->theta_vwp[2] = 0.20;    p->theta_vwp[3] = 0.20; p->theta_vwp[4] = 0.20;
		p->thermal_cond[0] = 5.8;  p->thermal_cond[1] = 5.8;   p->thermal_cond[2] = 5.7;      p->thermal_cond[3] = 5.8;  p->thermal_cond[4] = 5.8;
		p->psi_sat[0] = 0.26;   p->psi_sat[1] = 0.27;    p->psi_sat[2] = 0.28;      p->psi_sat[3] = 0.29;   p->psi_sat[4] = 0.31;
		break;

	case 8: // silty clay loam
		p->b[0] = 6.6;          p->b[1] = 6.8;           p->b[2] = 7.0;             p->b[3] = 7.2;          p->b[4] = 7.4;
		p->Ksat[0] = 0.00000042;    p->Ksat[1] = 0.00000038;     p->Ksat[2] = 0.00000034;       p->Ksat[3] = 0.000000252;   p->Ksat[4] = 0.000000084;
		p->fei[0] = 0.471;      p->fei[1] = 0.471;       p->fei[2] = 0.471;         p->fei[3] = 0.471;      p->fei[4] = 0.471;
		p->theta_vfc[0] = 0.37; p->theta_vfc[1] = 0.37;  p->theta_vfc[2] = 0.37;    p->theta_vfc[3] = 0.37; p->theta_vfc[4] = 0.37;
		p->theta_vwp[0] = 0.32; p->theta_vwp[1] = 0.32;  p->theta_vwp[2] = 0.32;    p->theta_vwp[3] = 0.32; p->theta_vwp[4] = 0.32;
		p->thermal_cond[0] = 4.2;  p->thermal_cond[1] = 4.2;   p->thermal_cond[2] = 4.2;     p->thermal_cond[3] = 4.2;  p->thermal_cond[4] = 4.2;
		p->psi_sat[0] = 0.33;   p->psi_sat[1] = 0.34;    p->psi_sat[2] = 0.35;      p->psi_sat[3] = 0.36;   p->psi_sat[4] = 0.38;
		break;


	case 9: // sandy clay
		p->b[0] = 6.0;          p->b[1] = 6.2;           p->b[2] = 6.4;             p->b[3] = 6.6;          p->b[4] = 6.8;
		p->Ksat[0] = 0.00000033;    p->Ksat[1] = 0.0000003;      p->Ksat[2] = 0.000000264;      p->Ksat[3] = 0.000000198;   p->Ksat[4] = 0.000000066;
		p->fei[0] = 0.430;      p->fei[1] = 0.430;       p->fei[2] = 0.430;         p->fei[3] = 0.430;      p->fei[4] = 0.430;
		p->theta_vfc[0] = 0.34; p->theta_vfc[1] = 0.34;  p->theta_vfc[2] = 0.34;    p->theta_vfc[3] = 0.34; p->theta_vfc[4] = 0.34;
		p->theta_vwp[0] = 0.24; p->theta_vwp[1] = 0.24;  p->theta_vwp[2] = 0.24;    p->theta_vwp[3] = 0.24; p->theta_vwp[4] = 0.24;
		p->thermal_cond[0] = 6.3;  p->thermal_cond[1] = 6.3;   p->thermal_cond[2] = 6.3;     p->thermal_cond[3] = 6.3;  p->thermal_cond[4] = 6.3;
		p->psi_sat[0] = 0.29;   p->psi_sat[1] = 0.30;    p->psi_sat[2] = 0.31;      p->psi_sat[3] = 0.32;   p->psi_sat[4] = 0.34;
		break;

	case 10: // silty clay
		p->b[0] = 7.9;          p->b[1] = 8.1;           p->b[2] = 8.3;             p->b[3] = 8.5;          p->b[4] = 8.7;
		p->Ksat[0] = 0.00000025;    p->Ksat[1] = 0.000000225;    p->Ksat[2] = 0.0000002;        p->Ksat[3] = 0.00000015;    p->Ksat[4] = 0.00000005;
		p->fei[0] = 0.479;      p->fei[1] = 0.479;       p->fei[2] = 0.479;         p->fei[3] = 0.479;      p->fei[4] = 0.479;
		p->theta_vfc[0] = 0.39; p->theta_vfc[1] = 0.39;  p->theta_vfc[2] = 0.39;    p->theta_vfc[3] = 0.39; p->theta_vfc[4] = 0.39;
		p->theta_vwp[0] = 0.25; p->theta_vwp[1] = 0.25;  p->theta_vwp[2] = 0.25;    p->theta_vwp[3] = 0.25; p->theta_vwp[4] = 0.25;
		p->thermal_cond[0] = 4.0;  p->thermal_cond[1] = 4.0;   p->thermal_cond[2] = 4.0;     p->thermal_cond[3] = 4.0;  p->thermal_cond[4] = 4.0;
		p->psi_sat[0] = 0.34;   p->psi_sat[1] = 0.35;    p->psi_sat[2] = 0.36;      p->psi_sat[3] = 0.37;   p->psi_sat[4] = 0.39;
		break;

	case 11: // clay
		p->b[0] = 7.6;          p->b[1] = 7.8;           p->b[2] = 8.0;             p->b[3] = 8.2;          p->b[4] = 8.4;
		p->Ksat[0] = 0.00000017;    p->Ksat[1] = 0.000000153;      p->Ksat[2] = 0.000000136;      p->Ksat[3] = 0.000000102;   p->Ksat[4] = 0.000000034;
		p->fei[0] = 0.475;      p->fei[1] = 0.475;       p->fei[2] = 0.475;         p->fei[3] = 0.475;      p->fei[4] = 0.475;
		p->theta_vfc[0] = 0.40; p->theta_vfc[1] = 0.40;  p->theta_vfc[2] = 0.40;    p->theta_vfc[3] = 0.40; p->theta_vfc[4] = 0.40;
		p->theta_vwp[0] = 0.27; p->theta_vwp[1] = 0.27;  p->theta_vwp[2] = 0.27;    p->theta_vwp[3] = 0.27; p->theta_vwp[4] = 0.27;
		p->thermal_cond[0] = 4.4;  p->thermal_cond[1] = 4.4;   p->thermal_cond[2] = 4.4;     p->thermal_cond[3] = 4.4;  p->thermal_cond[4] = 4.4;
		p->psi_sat[0] = 0.37;   p->psi_sat[1] = 0.38;    p->psi_sat[2] = 0.39;      p->psi_sat[3] = 0.40;   p->psi_sat[4] = 0.42;
		break;

	default:
		p->b[0] = 7.6;          p->b[1] = 7.8;           p->b[2] = 8.0;             p->b[3] = 8.2;          p->b[4] = 8.4;
		p->Ksat[0] = 0.00000017;    p->Ksat[1] = 0.000000153;      p->Ksat[2] = 0.000000136;      p->Ksat[3] = 0.000000102;   p->Ksat[4] = 0.000000034;
		p->fei[0] = 0.475;      p->fei[1] = 0.475;       p->fei[2] = 0.475;         p->fei[3] = 0.475;      p->fei[4] = 0.475;
		p->theta_vfc[0] = 0.40; p->theta_vfc[1] = 0.40;  p->theta_vfc[2] = 0.40;    p->theta_vfc[3] = 0.40; p->theta_vfc[4] = 0.40;
		p->theta_vwp[0] = 0.27; p->theta_vwp[1] = 0.27;  p->theta_vwp[2] = 0.27;    p->theta_vwp[3] = 0.27; p->theta_vwp[4] = 0.27;
		p->thermal_cond[0] = 4.4;  p->thermal_cond[1] = 4.4;   p->thermal_cond[2] = 4.4;     p->thermal_cond[3] = 4.4;  p->thermal_cond[4] = 4.4;
		p->psi_sat[0] = 0.37;   p->psi_sat[1] = 0.38;    p->psi_sat[2] = 0.39;      p->psi_sat[3] = 0.40;   p->psi_sat[4] = 0.42;

	}

}


// Init. the soil status: 
//	soil temperatures and moistures for each layer.
//	ponded water, snow depth et al.

// Liming He, Jun. 2, 2015

void Init_Soil_Status(struct Soil p[], double Tsoil, double Tair, double Ms, double snowdepth)
{
	int i;
	double d_t = Tsoil - Tair;

	p->Zp = 0.0; /* depth of ponded water on the surface*/
	p->Zsp = snowdepth;                /*Snow depth;*/
	p->r_rain_g = 0.0;  /*the rainfall rate, un--on understory g--on ground surface  m/s */


	if (d_t>5.0)  d_t = 5.0;
	if (d_t<-5.0) d_t = -5.0;

	p->temp_soil_c[0] = Tair + 0.4*d_t;
	p->temp_soil_c[1] = Tair + 0.5*d_t;
	p->temp_soil_c[2] = Tair + d_t;
	p->temp_soil_c[3] = Tair + 1.2*d_t;
	p->temp_soil_c[4] = Tair + 1.4*d_t;

	p->temp_soil_p[0] = Tair + 0.4*d_t;
	p->temp_soil_p[1] = Tair + 0.5*d_t;
	p->temp_soil_p[2] = Tair + d_t;
	p->temp_soil_p[3] = Tair + 1.2*d_t;
	p->temp_soil_p[4] = Tair + 1.4*d_t;

	//p->thetam[0] = 0.75*Ms;
	p->thetam[0] = 0.8*Ms;
	p->thetam[1] = Ms;
	p->thetam[2] = 1.05*Ms;
	p->thetam[3] = 1.10*Ms;
	p->thetam[4] = 1.15*Ms;

	//p->thetam_prev[0] = 0.75*Ms;
	p->thetam_prev[0] = 0.8*Ms;
	p->thetam_prev[1] = Ms;
	p->thetam_prev[2] = 1.05*Ms;
	p->thetam_prev[3] = 1.10*Ms;
	p->thetam_prev[4] = 1.15*Ms;

	for (i = 0; i < p->n_layer; i++)
	{
		if (p->temp_soil_c[i]<-1.0)  
			p->ice_ratio[i] = 1.0;
		else if (p->temp_soil_c[i]>0) 
			p->ice_ratio[i] = 0;
		else 
			p->ice_ratio[i] = (0 - p->temp_soil_c[i]) / 1.0;
	}
}

// funtion rewritten by LHE.  Jan 31, 2013.
// Last revision: May 16, 2015
// This function calculates the fraction of root in the soil for each soil layer. 

void SoilRootFraction(struct Soil soil[])
{
	int i;
	double cum_depth[MAX_LAYERS];

	// for the 0 layer
	cum_depth[0] = soil->d_soil[0];
	soil->f_root[0] = 1 - pow(soil->r_root_decay, cum_depth[0] * 100);

	// for 1 to n_layer-1
	for (i = 1; i < soil->n_layer - 1; i++)
	{
		cum_depth[i] = cum_depth[i - 1] + soil->d_soil[i];
		soil->f_root[i] = pow(soil->r_root_decay, cum_depth[i - 1] * 100) - pow(soil->r_root_decay, cum_depth[i] * 100);
	}

	// for the last layer. Put all reminding roots to the last layer. 
	// soil->Layer[layer].f_root = 1 - (1-pow(soil->r_root_decay, cum_depth[i-1] *100));
	soil->f_root[soil->n_layer - 1] = pow(soil->r_root_decay, cum_depth[soil->n_layer - 2] * 100);
}
