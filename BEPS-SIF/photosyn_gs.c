#include "beps.h"
#include "DB.h"

//juw
//Bin Chen: 28/05/2015
//Bo Qiu£º08/2019
/*

         This program solves a cubic equation to calculate
         leaf photosynthesis.  This cubic expression is derived from solving
         five simultaneous equations for A, PG, cs, CI and GS.
         Stomatal conductance is computed with the Ball-Berry model.
         The cubic derivation assumes that b', the intercept of the Ball-Berry
         stomatal conductance model, is non-zero.

          Gs = k A rh/cs + b'


          We also found that the solution for A can be obtained by a quadratic equation
          when Gs is constant or b' is zero.


            The derivation is published in:

            Baldocchi, D.D. 1994. An analytical solution for coupled leaf photosynthesis
            and stomatal conductance models. Tree Physiology 14: 1069-1079.


-----------------------------------------------------------------------

          A Biochemical Model of C3 Photosynthesis

            After Farquhar, von Caemmerer and Berry (1980) Planta.
            149: 78-90.

        The original program was modified to incorporate functions and parameters
        derived from gas exchange experiments of Harley, who paramertized Vc and J in
        terms of optimal temperature, rather than some reference temperature, eg 25C.

        Program calculates leaf photosynthesis from biochemical parameters

        rd25 - Dark respiration at 25 degrees C (umol m-2 s-1)
        tlk - leaf temperature, Kelvin
        jmax - optimal rate of electron transport
        vcopt - maximum rate of RuBP Carboxylase/oxygenase
        iphoton - incident photosynthetically active photon flux (umols m-2 s-1)
   
            note: Harley parameterized the model on the basis of incident PAR

        gs - stomatal conductance (mols m-2 s-1), typically 0.01-0.20
        pstat-station pressure, bars
        aphoto - net photosynthesis  (umol m-2 s-1)
        ps - gross photosynthesis (umol m-2 s-1)
        aps - net photosynthesis (mg m-2 s-1)
        aphoto (umol m-2 s-1)

--------------------------------------------------

        iphoton is radiation incident on leaves

        The temperature dependency of the kinetic properties of
        RUBISCO are compensated for using the Arrhenius and
        Boltzmann equations.  From biochemistry, one observes that
        at moderate temperatures enzyme kinetic rates increase
        with temperature.  At extreme temperatures enzyme
        denaturization occurs and rates must decrease.

        Arrhenius Eq.

        f(T)=f(tk_25) exp(tk -298)eact/(298 R tk)), where eact is the
        activation energy.

        Boltzmann distribution

        F(T)=tboltz)


        Define terms for calculation of gross photosynthesis, PG

        PG is a function of the minimum of RuBP saturated rate of
        carboxylation, Wc, and the RuBP limited rate of carboxylation, Wj.
        Wj is limiting when light is low and electron transport, which
        re-generates RuBP, is limiting.  Wc is limiting when plenty of RuBP is
        available compared to the CO2 that is needed for carboxylation.

        Both equations take the form:

        PG-photorespiration= (a CI-a d)/(e CI + b)

        PG-photorespiration=min[Wj,Wc] (1-gamma_ps/Ci)

        Wc=Vcmax Ci/(Ci + Kc(1+O2/Ko))

        Wj=J Ci/(4 Ci + 8 gamma_ps)

        Ps kinetic coefficients from Harley at WBW.

        gamma_ps is the CO2 compensation point

        Jan 14, 1999 Updated the cubic solutions for photosynthesis.  There are
        times when the restriction that R^2 < Q^3 is violated.  I therefore need
        alternative algorithms to solve for the correct root.

        Solar-induced chlorophyll fluorescence at the leaf level was simulated. 
==============================================================*/

void photosynthesis(lai_leaf_p,temp_leaf_p,rad_leaf,e_air,g_lb_w,vc_opt,f_soilwater,b_h2o,m_h2o,cii,temp_leaf_c,LH_leaf,Gs_w,aphoto,ci,xSIF)

double lai_leaf_p;
double temp_leaf_p; // It is temporay, will be removed later

double rad_leaf; // net shortwave radiation (W/m2)
double e_air;  // water vapor pressure above canopy (kPa)
double g_lb_w; // leaf laminar boundary layer condunctance to H2O (m/s)
double vc_opt; // the maximum velocities of carboxylation of Rubisco at 25 deg C (umol m-2 s-1)
		
double f_soilwater; // an empirical scalar of soil water stress on stomatal conductance, dimensionless
double b_h2o; // the intercept term in BWB model (mol H2O m-2 s-1)
double m_h2o; // the slope in BWB model


double cii; // initial intercellular co2 concentration (ppm)
double temp_leaf_c; // leaf temperature (deg C)

double LH_leaf; // leaf latent heat flux (W m-2)

double *Gs_w; // stomatal conductance to water vapor (m s-1)
double *aphoto; // net photosynthesis rate (umol CO2 m-2 s-1)
double *ci; // intercellular co2 concentration (ppm)
//QBO
double *xSIF;

{
	double air_pres=101.325; // air pressure (kPa)
	double ca;  // atmospheric co2 concentration (ppm)
	double iphoton; // incident photosynthetic photon flux density (PPFD) umol m-2 s-1
	double g_lb_c; // leaf laminar boundary layer condunctance to CO2 (mol m-2 s-1)
	double rh_leaf; // relative humidity at leaf surface (0-1)
	double temp_leaf_K; // leaf temperature (K)
	double gs_co2_mole; // stomatal conductance to CO2 (mol m-2 s-1)
	double gs_h2o_mole; // stomatal conductance to h2o (mol m-2 s-1)
        double bc; // temporary variable
	double cs; // CO2 concentration at leaf surface (ppm)		

		double b_co2; // the intercept term in BWB model (mol CO2 m-2 s-1): b_h2o/1.6
		double m_co2; // the slope in BWB model: m_h2o/1.6
	
		double gammac; //CO2 compensation point (ppm)
		double jmopt; //the maximum potential electron transport rate at 25 deg C (umol m-2 s-1)
        double jmax; //the maximum potential electron transport rate (umol m-2 s-1)
		double vcmax; //the maximum velocities of carboxylation of Rubisco (umol m-2 s-1)
        double km_co2; // Michaelis-Menten constant for CO2 (µmol mol-1) 
		double km_o2; // Michaelis-Menten constant for O2 (mmol mol-1)
		double tau; // the specifity of Rubisco for CO2 compared with O2 
        double resp_ld; // leaf dark respiration (umol m-2 s-1)
		double resp_ld25; // leaf dark respiration at 25 deg C (umol m-2 s-1)
		
		double j_photon; //the flux of electrons through the thylakoid membrane (umol m-2 s-1) 
		double alpha_ps;
		double beta_ps;
		double gamma_ps;
		double theta_ps;

        double denom;
		double p_cubic;
		double q_cubic;
		double r_cubic;
        
		double Qroot;
		double Rroot;
        double root1,root2,root3;
        double ang_L;
       
		double j_sucrose; // net photosynthesis rate limited by sucrose synthesis (umol m-2 s-1)
		double wc,wj,psguess; // gross photosynthesis rate limited by light (umol m-2 s-1)

        double je,xxn,ffs;   //QBO
	
        double Aquad,Bquad,Cquad;
        double b_ps, a_ps, e_ps, d_ps;		
        double product;
		double ps_1;
		double delta_1;
        double r3q;
        double minroot=0, maxroot=0, midroot=0;
		
		double tprime25;

	ca=CO2_air;	
		iphoton = 4.55*0.5*rad_leaf; // 

		if(2*iphoton < 1) 
			iphoton = 0;

		temp_leaf_K = temp_leaf_c + 273.13;

		fact.latent		= LAMBDA(temp_leaf_p);
		bound_layer_res.vapor = 1.0/g_lb_w;

//		g_lb_c = (g_lb_w/1.6)*air_pres/(temp_leaf_K*rugc); // (mol m-2 s-1)

		met.press_bars	= 1.013;
		met.pstat273	= 0.022624 / (273.16 * met.press_bars);

		met.T_Kelvin		= temp_leaf_c+273.13;
		met.rhova_g		= e_air * 2165/met.T_Kelvin;	// absolute humidity, g m-3  
		met.rhova_kg		= met.rhova_g / 1000.;		// absolute humidity, kg m-3
		
		g_lb_c	= 1. / (1.0/g_lb_w*1.6 * temp_leaf_K * (met.pstat273));

        m_co2 = m_h2o/1.6;
        b_co2 = b_h2o/1.6;

		rh_leaf  = SFC_VPD(temp_leaf_K, LH_leaf);	
		
		tprime25	= temp_leaf_K - tk_25;     // temperature difference

        /** use Arrhenius Eq. to compute KC and km_o2 */
	km_co2	= TEMP_FUNC(kc25, ekc, tprime25, tk_25, temp_leaf_K);
	km_o2	= TEMP_FUNC(ko25, eko, tprime25, tk_25,temp_leaf_K);
	tau	= TEMP_FUNC(tau25, ektau, tprime25, tk_25,temp_leaf_K);

        bc = km_co2 * (1.0 + o2 / km_o2);

  /*      if(iphoton < 1)
			iphoton = 0;*/

        gammac = 0.5 * o2/tau*1000.0; // umol mol-1

        resp_ld25=vc_opt * 0.004657;

        if(2.0*iphoton > 10) //Bin Chen: check this later. reduce respiration by 40% in light according to Amthor 
			resp_ld25 *= 0.4;      

		resp_ld = TEMP_FUNC(resp_ld25, erd, tprime25, tk_25, temp_leaf_K);
        
//		jmopt = 29.1 + 1.64*vc_opt;
	
		jmopt		= 2.39*vc_opt - 14.2; 


        jmax = TBOLTZ(jmopt, ejm, toptjm, temp_leaf_K); // Apply temperature correction to JMAX     
        vcmax = TBOLTZ(vc_opt, evc, toptvc, temp_leaf_K);// Apply temperature correction to vcmax


/***************************************

        APHOTO = PG - resp_ld, net photosynthesis is the difference
        between gross photosynthesis and dark respiration. Note
        photorespiration is already factored into PG.

****************************************

/*
        Gs from Ball-Berry is for water vapor.  It must be divided
        by the ratio of the molecular diffusivities to be valid
        for A
*/
        
        alpha_ps = 1.0 + (b_co2 / g_lb_c) - m_co2*rh_leaf*f_soilwater;
        beta_ps = ca * (g_lb_c*m_co2*rh_leaf*f_soilwater - 2.0 * b_co2 - g_lb_c);
        gamma_ps = ca * ca * g_lb_c * b_co2;
        theta_ps = g_lb_c*m_co2*rh_leaf*f_soilwater - b_co2;

/*
        Test for the minimum of Wc and Wj.  Both have the form:

        W = (a ci - ad)/(e ci + b)

        after the minimum is chosen set a, b, e and d for the cubic solution.

        estimate of J according to Farquhar and von Cammerer (1981)
*/
            
 /*     if (jmax > 0)
			j_photon = qalpha * iphoton / sqrt(1. +(qalpha2 * iphoton * iphoton / (jmax * jmax)));
        else
			j_photon = 0;*/  //J photon from Harley

		j_photon =jmax * iphoton/ (iphoton+ 2.1*jmax);		


	// initial guess of intercellular CO2 concentration to estimate Wc and Wj:
        wj = j_photon * (cii - gammac) / (4. * cii + 8.0*gammac);

        wc = vcmax * (cii - gammac) / (cii + bc);

        if(wj < wc)
        {

        // for Harley and Farquhar type model for Wj

			psguess=wj;
			a_ps = j_photon;        
			b_ps = 8.0*gammac;        
			e_ps = 4.0;
			d_ps = gammac;
        }
        else
        {
			psguess=wc;
			a_ps = vcmax;
			b_ps = bc;			
			e_ps = 1.0;
			d_ps = gammac;
        }

        //if wj or wc are less than resp_ld then A would probably be less than zero.  This would yield a
       // negative stomatal conductance.  In this case, assume gs equals the cuticular value. This
       // assumptions yields a quadratic rather than cubic solution for A
 

        if (wj <= resp_ld)
			goto quad;

        if (wc <= resp_ld)
			goto quad;

        //cubic solution:
        // A^3 + p A^2 + q A + r = 0

        denom = e_ps * alpha_ps;

        p_cubic = (e_ps * beta_ps + b_ps * theta_ps - a_ps * alpha_ps + e_ps * resp_ld * alpha_ps);
        p_cubic /= denom;

        q_cubic = (e_ps * gamma_ps + (b_ps * gamma_ps / ca) - a_ps * beta_ps + a_ps * d_ps * theta_ps + e_ps * resp_ld * beta_ps + resp_ld * b_ps * theta_ps);
        q_cubic /= denom;

        r_cubic = -a_ps * gamma_ps + a_ps * d_ps * gamma_ps/ca + e_ps * resp_ld * gamma_ps + resp_ld * b_ps * gamma_ps/ca;
        r_cubic /= denom;


        // Use solution from Numerical Recipes from Press

        
        Qroot = (p_cubic*p_cubic - 3.0 * q_cubic) / 9.0;
        Rroot = (2.0 * p_cubic*p_cubic*p_cubic  - 9.0 * p_cubic * q_cubic + 27.0 * r_cubic) / 54.0;
  
       	r3q = Rroot / sqrt(Qroot*Qroot*Qroot);
	if (r3q>1) r3q=1;	//  by G. Mo
	if (r3q<-1) r3q=-1;	//  by G. Mo

	ang_L = acos(r3q);
	
        root1 = -2.0 * sqrt(Qroot) * cos(ang_L / 3.0) - p_cubic / 3.0;  // real roots 
        root2 = -2.0 * sqrt(Qroot) * cos((ang_L + PI2) / 3.0) - p_cubic / 3.0;
        root3 = -2.0 * sqrt(Qroot) * cos((ang_L - PI2) / 3.0) - p_cubic / 3.0;

		
        // Here A = x - p / 3, allowing the cubic expression to be expressed
         //as: x^3 + ax + b = 0

    // rank roots #1,#2 and #3 according to the minimum, intermediate and maximum value

	   if(root1 <= root2 && root1 <= root3)
       {
		   minroot=root1;           
		   if (root2 <= root3)
           {
			   midroot=root2;
               maxroot=root3;
		   }
           else
           {
			   midroot=root3;
               maxroot=root2;
		   }
	   }


       if(root2 <= root1 && root2 <= root3)
       {
		   minroot=root2;           
		   if (root1 <= root3)
           {
			   midroot=root1;
               maxroot=root3;
		   }
           else
		   {
			   midroot=root3;
               maxroot=root1;
		   }
       }


       if(root3 <= root1 && root3 <= root2)
       {
		   minroot=root3;
           if (root1 < root2)
           {
			   midroot=root1;
               maxroot=root2;
		   }
           else
           {
			   midroot=root2;
               maxroot=root1;
		   }

        }  

		*aphoto=0;
        // find out where roots plop down relative to the x-y axis
            
        
        if (minroot > 0 && midroot > 0 && maxroot > 0)
			*aphoto=minroot;
        
        
        if (minroot < 0 && midroot < 0 && maxroot > 0)
			*aphoto=maxroot;

   
        if (minroot < 0 && midroot > 0 && maxroot > 0)
			*aphoto=midroot;

 
        //also test for sucrose limitation of photosynthesis, as suggested by
        //Collatz.  Js=Vmax/2
 
        j_sucrose = vcmax / 2. - resp_ld;

        if(j_sucrose < *aphoto)
			*aphoto = j_sucrose;

       
         //Stomatal conductance for water vapor

        //forest are hypostomatous.
        //Hence we don't divide the total resistance
        //by 2 since transfer is going on only one side of a leaf.
  
       
        // if A < 0 then gs should go to cuticular value and recalculate A
         //using quadratic solution
        
        if(*aphoto <= 0.0)
			goto quad;
		else
			goto OUTDAT;



        // if aphoto < 0  set stomatal conductance to cuticle value   

quad:
        //a quadratic solution of A is derived if gs=b, but a cubic form occurs
        //if gs =ax + b.  Use quadratic case when A <=0 
         

  /* Bin Chen:
         r_tot = 1.0/b_co2 + 1.0/g_lb_c; // total resistance to CO2 (m2 s mol-1)
         denom = g_lb_c * b_co2;

         Aquad = r_tot * e_ps;
         Bquad = (e_ps*resp_ld + a_ps)*r_tot - b_ps - e_ps*ca;
         Cquad = a_ps*(ca-d_ps) - resp_ld*(e_ps*ca+b_ps);

*/

		// original version

		ps_1	= ca * g_lb_c * b_co2;
		delta_1	= b_co2 + g_lb_c;
		denom	= g_lb_c * b_co2;
	
		Aquad = delta_1 * e_ps;
		Bquad = -ps_1 * e_ps - a_ps * delta_1 + e_ps * resp_ld * delta_1 - b_ps * denom;
		Cquad = a_ps * ps_1 - a_ps * d_ps * denom - e_ps * resp_ld * ps_1 - resp_ld * b_ps * denom;


        product=Bquad * Bquad - 4.0 * Aquad * Cquad;

        if (product >= 0)
		//	*aphoto = (-Bquad + sqrt(product)) / (2.0 * Aquad);
			*aphoto = (-Bquad - sqrt(product)) / (2.0 * Aquad);   
       
OUTDAT:
		*aphoto =	max(0, *aphoto);

        cs = ca - *aphoto / g_lb_c;       

		gs_h2o_mole = (f_soilwater *m_h2o* rh_leaf * *aphoto / cs) + b_h2o; // mol m-2 s-1
		gs_co2_mole = gs_h2o_mole /1.6;

		*ci = cs - *aphoto / gs_co2_mole;

		*Gs_w = gs_h2o_mole * temp_leaf_K * (met.pstat273); // m s-1   

//QBO

        je=*aphoto*(*ci+2.0*gammac)/(*ci-gammac);
        if (je<0)
        je = 0.0;

        xxn = 1.0-je/(iphoton*0.05);

        if (xxn < 0)
        xxn = 0.0;

        ffs = SIF_y(xxn);
        *xSIF = ffs*iphoton;
 
        return;
}

/*----------------------------------------------------------------*/
 double SIF_y (double xxn)
/*----------------------------------------------------------------*/
/* this function computes the solar-induced fluorescence 
   at the leaf level, by Bo Qiu, 08/2015 
 *  *
 *   *       ----------------------------------------------------------------*/
 {
         double y, kf, kd, kp, kk;
         double ps_sif, kn, fm, fs;

         kf = 0.05;
         kd = 0.95;
         kp = 4.0;

         ps_sif = kp/(kf+kp+kd)*(1.0-xxn);
//         kn = (6.2473*xxn-0.5994)*xxn;
         kk=exp(log(xxn)*2.83);
         kn = 2.48*(1+0.114)*kk/(kk+0.114);
         fm = kf/(kf+kd+kn);
         fs = fm*(1.0-ps_sif);
         y  = fs;

         return y;
 }

/*----------------------------------------------------------------*/
 double SFC_VPD (double temp_leaf_K, double leleafpt)

/* this function computes the relative humidity at the leaf surface for
 application in the Ball Berry Equation
 latent heat flux, LE, is passed through the function, mol m-2 s-1
 and it solves for the humidity at leaf surface

 es_leaf : saturation vapor pressure at leaf temperature.

  ----------------------------------------------------------------*/
 {
	 double y, rhov_sfc,e_sfc,vpd_sfc,rhum_leaf;
	 double es_leaf;
	 
	 es_leaf  = ES(temp_leaf_K);
	 rhov_sfc = (leleafpt / (fact.latent)) * bound_layer_res.vapor + met.rhova_kg;  /* kg m-3 */
	 
	 e_sfc = rhov_sfc * temp_leaf_K / .2165;			// mb 
	 vpd_sfc = es_leaf - e_sfc;			// mb 
	 rhum_leaf = 1. - vpd_sfc / es_leaf;		// 0 to 1.0 
	 y = rhum_leaf;
	 
	 return y;
 }

/*-------------------------------------------------------------------------------*/
 double TEMP_FUNC(double rate,double eact,double tprime,double tref, double t_lk)
/*-------------------------------------------------------------------------------*/
 {

        //  Arhennius temperature function  
            
        double y;
        y = rate * exp(tprime * eact / (tref * rugc*t_lk));
        return y;
 }

/*----------------------------------------------------------------*/
 double LAMBDA (double tak)
/*----------------------------------------------------------------*/
 {
	  // Latent heat of Vaporiation, J kg-1
	  
	  double y;
	  
	  y = 3149000. - 2370. * tak;
	  
	  // add heat of fusion for melting ice  
	  if(tak < 273.)
		  y += 333;
	  
	  return y;
 }

/*----------------------------------------------------------------*/
 double ES(double t)

// saturation vapor pressure function (mb)
// T is temperature in Kelvin
/*----------------------------------------------------------------*/
 {
	 double y, y1;
	 
	 if(t > 0)       
	 {
		 y1 = (54.8781919 - 6790.4985 / t - 5.02808 * log(t));
		 y = exp(y1);
	 }
	 else
		 printf("bad es calc");
	 
	 return y;
 }

/*----------------------------------------------------------------*/
 double TBOLTZ(double rate, double eakin, double topt, double tl)
// Boltzmann temperature distribution for photosynthesis 
/*----------------------------------------------------------------*/
 {
	 
	 double y, dtlopt,prodt,numm,denom;
	 
	 dtlopt	= tl   - topt;
	 prodt	= rugc * topt * tl;
	 numm	= hkin * exp(eakin * (dtlopt) / (prodt));
	 denom	= hkin - eakin * (1.0 - exp(hkin * (dtlopt) / (prodt)));
	 y	= rate * numm / denom;
	 return y;
 }

