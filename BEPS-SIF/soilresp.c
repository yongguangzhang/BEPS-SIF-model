/*************************************************************************
    soilresp.c

	Calculating soil respiration.
*****************************************************************************/
#include "beps.h"
#include "soil.h"
 

void soilresp(Ccd,Cssd,Csmd,Cfsd,Cfmd,Csm,Cm,Cs,Cp,npp_yr,coef,soiltype,soilp,mid_res)

double *Ccd;
double *Cssd;
double *Csmd;
double *Cfsd;
double *Cfmd;
double *Csm;
double *Cm;
double *Cs;
double *Cp;
float npp_yr;
double *coef;
int soiltype;
struct Soil *soilp;
struct results *mid_res;


{  
	double fw, fcr, fl, ffr; 
	double kw_cd, kcr_cd;
	double kl_sl, kfr_fl; 
	double km_p, ks_p;
   	double kssd_a, kssd_sm, kssd_s, ksmd_a, ksmd_sm,kfsd_a, kfsd_m, kfsd_s, kfmd_a, kfmd_m;
	double kcd_a, kcd_m;
   	double kcd_s,ksm_a,ksm_s, km_a, km_s, ks_a, ks_m,kp_a, kp_m; 	
	double Cw[10], Ccr[10];
	double Cl[10], Cfr[10];
	double dCw[10], dCcr[10];
	double dCl[10], dCfr[10];
	double dCcd[10],dCssd[10],dCsmd[10],dCfsd[10],dCfmd[10];
	double dCsm[10], dCm[10];
	double dCs[10], dCp[10];
	double part1,part2;	
	double Fm[10],npp;
	double lambda[layer+1],lambda_t[layer+1],lambda_w[layer+1];	
	double lam_u,lam_d;
	int  ii;

	for(ii=1;ii<=layer;ii++)
	{
		 lambda_t[ii]=exp(308.56*(1/(35.0+46.032)-1/(46.032+soilp->temp_soil_c[ii - 1])));
		 lambda_t[ii]=min(1.0,lambda_t[ii]);
		 lambda_t[ii]=max(0.3,lambda_t[ii]);
	}
	for(ii=1;ii<=layer;ii++)
	{
		 if(soiltype>=6)	
			lambda_w[ii]=5.44*soilp->thetam[ii-1]/soilp->fei[ii-1]-5.03*pow(soilp->thetam[ii-1]/soilp->fei[ii-1],2)-0.472; 
		 else 
			lambda_w[ii]=5.63*soilp->thetam[ii-1]/soilp->fei[ii-1]-4.64*pow(soilp->thetam[ii-1]/soilp->fei[ii-1],2)-0.710; 
		 lambda_w[ii]=max(0.3,lambda_w[ii]); 
	}
	for(ii=1;ii<=layer;ii++)  lambda[ii]=lambda_t[ii]*lambda_w[ii];
	
	lam_u= lambda[1];                    /*for surface pool*/
    	lam_d= lambda[2];    /* for soil pool */
	
    	fw = coef[0];
	fcr = coef[1];
	fl = coef[2];
	ffr = coef[3];
	kw_cd =coef[4]/8760 ;
	kcr_cd =coef[5]/8760;
	kl_sl =coef[6]/8760;
	kfr_fl =coef[7]/8760; 
	kssd_a = coef[8]/8760;
   	kssd_sm = coef[9]/8760;      
    	kssd_s = coef[10]/8760;
	ksmd_a = coef[11]/8760;
	ksmd_sm = coef[12]/8760;
    	kfsd_a = coef[13]/8760;
    	kfsd_m = coef[14]/8760;      
    	kfsd_s = coef[15]/8760;
	kfmd_a = coef[16]/8760;
	kfmd_m = coef[17]/8760;	
	kcd_a = coef[18]/8760;
    	kcd_m = coef[19]/8760;
    	kcd_s = coef[20]/8760;
    	km_a = coef[21]/8760;
    	km_p = coef[22]/8760;
	km_s = coef[23]/8760;
    	ksm_a= coef[24]/8760;
	ksm_s= coef[25]/8760;
	ks_a = coef[26]/8760;
	ks_p = coef[27]/8760;
    	ks_m = coef[28]/8760;
	kp_a = coef[29]/8760;
    	kp_m = coef[30]/8760;

	
	Cw[0]   =  coef[0]/coef[4] * npp_yr;  // for stem  gC/m2
	Ccr[0]  =  coef[1]/coef[5] * npp_yr;  // for coast root   gC/m2
	Cl[0]   =  coef[2]/coef[6] * npp_yr;  // for leaf   gC/m2
	Cfr[0]  =  coef[3]/coef[7] * npp_yr;  // for fine root   gC/m2
	
	Fm[1]=0.2;

	npp = mid_res->npp_o + mid_res->npp_u;
	
	dCw[1]  = fw * npp- kw_cd  * Cw[0];	
	dCcr[1] = fcr* npp - kcr_cd * Ccr[0]; 
	dCl[1]  = fl * npp - kl_sl  * Cl[0];
    	dCfr[1] = ffr* npp- kfr_fl * Cfr[0]; 
	
	Cw[1]  =Cw[0] + dCw[1];
	Ccr[1] =Ccr[0] + dCcr[1];
	Cl[1]  =Cl[0] + dCl[1];
	Cfr[1] =Cfr[0] + dCfr[1];  

	part1=(kw_cd * Cw[1]+kcr_cd * Ccr[1])/(1+lam_d*(kcd_a + kcd_m + kcd_s));
	part2=Ccd[0] * lam_d* (kcd_a + kcd_m + kcd_s);
	dCcd[1]=part1-part2;
	Ccd[1] = Ccd[0] + dCcd[1];   
	/* Coarse detritus from woody and coarse root;*/
	
	part1=(1 - Fm[1])* kl_sl*Cl[1]/(1+lam_u*(kssd_a + kssd_sm + kssd_s));
	part2=Cssd[0]* lam_u * (kssd_a + kssd_sm + kssd_s);
	dCssd[1]=part1-part2;
	Cssd[1]= Cssd[0]+dCssd[1];      
	/* for surface structural litter*/
	
	part1=Fm[1]* kl_sl * Cl[1]/(1+lam_u*(ksmd_a + ksmd_sm));
	part2= Csmd[0]* lam_u * (ksmd_a + ksmd_sm);
	dCsmd[1]=part1-part2;
	Csmd[1]= Csmd[0]+dCsmd[1];         
	/* for surface metabolic litter*/
	
	part1=(1 - Fm[1])* kfr_fl* Cfr[1]/(1+lam_d*(kfsd_a + kfsd_m + kfsd_s));
	part2=Cfsd[0]* lam_d * (kfsd_a + kfsd_m + kfsd_s);
	dCfsd[1]=part1-part2;
	Cfsd[1]= Cfsd[0]+dCfsd[1];          
	/*for soil strutural litter pool*/
	
	part1=Fm[1] * kfr_fl * Cfr[1]/(1+lam_d * (kfmd_a + kfmd_m)); 
	part2=lam_d * (kfmd_a + kfmd_m)* Cfmd[0];
	dCfmd[1]=part1-part2;
	Cfmd[1]= Cfmd[0]+dCfmd[1];       
	/* for soil metabolic pool*/
	
	part1=lam_u*(Cssd[1]*kssd_sm+Csmd[1]*ksmd_sm);
	part2=lam_u*Csm[0]*(ksm_a+ksm_s);
	dCsm[1]=part1-part2;
	Csm[1]=Csm[0]+dCsm[1];                     
	/* for surface microbe pool*/
	
	part1=(lam_d * (kfsd_m * Cfsd[1]+kfmd_m*Cfmd[1] + Ccd[1] * kcd_m) +lam_d*(Cs[1-1]*ks_m+Cp[1-1] * kp_m));				
	part2=Cm[0] * lam_d*(km_a +  km_s +km_p);
	dCm[1]=part1-part2;
	Cm[1]=Cm[0]+dCm[1];                          
	/* for soil microbe pool*/
	   
	part1=(lam_d*(Cm[1]*km_s + Ccd[1] * kcd_s +Cfsd[1]*kfsd_s )+ lam_u* (Csm[1]*ksm_s + Cssd[1]*kssd_s));
	part2=Cs[0]* lam_d *( ks_a + ks_p+ks_m);
	dCs[1]=part1-part2;
	Cs[1]=Cs[0]+dCs[1];   
	/* for slow carbon pool*/
	
	dCp[1] =(lam_d *( km_p * Cm[1] + ks_p * Cs[1]) - lam_d * (kp_m * Cp[1-1] + kp_a * Cp[1-1]));
	Cp[1]=Cp[0]+dCp[1];                             
	/* for passive carbon pool.*/
	
	
	/** NEP  ***/
	mid_res->NEP = npp+(dCsmd[1]+dCssd[1]+dCfsd[1]+dCfmd[1]+dCcd[1]+dCm[1]+dCsm[1]+dCs[1]+dCp[1]);
	/** NPP  ***/
	mid_res->NPP= npp;
	
	/* soil respiration */
//	mid_res->soil_resp = npp - mid_res->NEP;

    return;
} 
