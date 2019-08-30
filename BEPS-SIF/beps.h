/*************************************************************************
  program: 	beps.h
  Description:  Header file for defining constants and global variables
  ----------    for BEPS program

***************************************************************************
  CCRS (EMS/Applications Division)
  Wrintten by: 	J. Liu
  Modified by:  G. Mo
  Last update:  June 2015
*****************************************************************************/

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<malloc.h>
#include <string.h>

/* Define Constants */
#define NOERROR		0
#define ERROR		1
#define PI			3.1415926
#define zero		0.0000000001

//#define max(a,b) 	(a>b)?a:b       // used for UNIX
//#define min(a,b)	(a<b)?a:b	// used for UNIX
#define max(a,b) 	((a)>(b))?(a):(b) // LHE. the  orginal one can lead to disorder. 
#define min(a,b)	((a)<(b))?(a):(b) // LHE


#define l_sta		105  // start line
#define l_end		105  // end line
#define p_sta		101  // start pix
#define p_end		101  // end pix

#define RTIMES 		24      /*24  */
#define step 	3600	/*3600 in sec */
#define kstep		371.14        /* 10 times per hour, 360 sec. per time */
#define kloop		10        /* 10 times per hour, 360 sec. per time */
#define layer		5
#define depth_f		6
#define  CO2_air         369.54      //atmopsheric CO2 concentration
#define rho_a	1.292   // densoty of air at 0C


void readconf();
void mid_prg();
void readinput1();
void readinput2();
void readlai_d();
void readlonlat();
void readlai_all();
void inter_prg();
void s_coszs();
void aerodynamic_conductance();
void plantresp();
void Vcmax_Jmax();
void netRadiation ();
void soilresp();
void readparam();
void lai2();
void readcoef();
void photosynthesis();
void Leaf_Temperatures();
double Leaf_Temperature();

void sensible_heat();
void transpiration();
void evaporation_canopy();
void evaporation_soil();
void rainfall_stage1();
void rainfall_stage2();
void rainfall_stage3();

void meteo_pack();
void surface_temperature();
void snowpack_stage1();
void snowpack_stage2();
void snowpack_stage3();

/*	Declare structures */
struct climatedata
	{
	double Srad;
	double LR;
	double temp;
	double rh;
	double rain;
	double wind;
	double dr_o;
	double df_o;
	double dr_u;
	double df_u;
//	float st_c;
	};

struct results
	{
	double gpp_o_sunlit;
	double gpp_u_sunlit;
	double gpp_o_shaded;
	double gpp_u_shaded;
	double plant_resp;
	double npp_o;
	double npp_u;
	double GPP;
	double NPP;
	double NEP;
	double soil_resp;
	double Net_Rad;
	double SH;
	double LH;
	double Trans;
	double Evap;
        double SIF;
	};


/*	Declare global variables */
short lc_no;
int yr,bgn_day,end_day;
int npixels,nlines;

char lc_fn[255];   /* Land cover file */
char lai_fn[255];    /* Leaf area index file  */
char lai_fp[255];    /* Leaf area index file prefix */
char stxt_fn[255];         /* soil texture file */
char ci_fn[255];    /* clumping index file */
char st_fn[255];         /* intial values of soil temp */
char sw_fn[255];         /* intial values of soil water */
char sdp_fn[255];         /* intial values of snow depth*/

char r_fn[255];         /* meteo.data files */
char t_fn[255];
char h_fn[255];
char p_fn[255];
char wd_fn[255];

char lon_fn[255];
char lat_fn[255];

char laiy_fn[255],nppy_fn[255];
char ccd_fn[255],cssd_fn[255],csmd_fn[255],cfsd_fn[255],cfmd_fn[255],csm_fn[255],cm_fn[255],cs_fn[255],cp_fn[255];

char fp4outp1[255];    /* output file1 prefix */
char fp4outp2[255];    /* output file2 prefix */
char fp4outp3[255];    /* output file3 prefix */

