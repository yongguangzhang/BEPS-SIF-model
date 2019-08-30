
#define PI180 0.017453292                       // pi divided by 180, radians per degree
#define PI9 2.864788976
#define PI2 6.283185307                         // 2 time pi

 struct meteorology {

                double ustar;                   // friction velocity, m s-1
                double ustarnew;                // updated friction velocity with new H, m s-1
                double rhova_g;                 // absolute humidity, g m-3
                double rhova_kg;                // absolute humidity, kg m-3
                double sensible_heat_flux;      // sensible heat flux, W M-2
                double H_old;                   // old sensible heat flux, W m-2              
                double air_density;             // air density, kg m-3
                double T_Kelvin;                // absolute air temperature, K
                double press_kpa;               // station pressure, kPa
                double press_bars;              // station pressure, bars
                double press_Pa;                // pressure, Pa
                double pstat273;                // gas constant computations
                double air_density_mole;        // air density, mole m-3
                double relative_humidity;       // relative humidity, ea/es(T)
                double vpd;                     // vapor pressure deficit 
				double ir_in;                   // infrared flux density
                } met;

 // structure for plant and physical factors

                struct factors {
                double latent;      // latent heat of vaporization, J kg-1
                double latent18;    // latent heat of vaporization times molecular mass of vapor, 18 g mol-1
                double heatcoef;    // factor for sensible heat flux density
                double a_filt;          // filter coefficients
                double b_filt;      // filter coefficients
                double co2;         // CO2 factor, ma/mc * rhoa (mole m-3)
        
        } fact;

 struct boundary_layer_resistances{
                
                double vapor;                   // resistance for water vapor, s/m
                double heat;                    // resistance for heat, s/m
                double co2;                     // resistance for CO2, s/m
                } bound_layer_res;

void   TBOLTZdouble();
double TEMP_FUNC();
double TBOLTZ();
void   photosynthesis();
double SFC_VPD();
double ES();
double LAMBDA();
double SIF_y();




#define  rugc   8.314              // J mole-1 K-1 
//#define  vcopt  73.0   // carboxylation rate at optimal temperature, umol m-2 s-1 
//#define  jmopt  170.0  // electron transport rate at optimal temperature, umol m-2 s-1 
#define  rd25   0.34     // dark respiration at 25 C, rd25= 0.34 umol m-2 s-1 
#define  pi4    12.5663706
        //  Universal gas constant  
        
#define  rgc1000 8314            // gas constant times 1000.

        // Consts for Photosynthesis model and kinetic equations.
        // for Vcmax and Jmax.  Taken from Harley and Baldocchi (1995, PCE)
#define hkin  200000.0   // enthalpy term, J mol-1
#define skin  710.0      // entropy term, J K-1 mol-1
#define ejm   55000.0      // activation energy for electron transport, J mol-1
#define evc   55000.0      // activation energy for carboxylation, J mol-1
    
        // Enzyme constants & partial pressure of O2 and CO2
        // Michaelis-Menten K values. From survey of literature.

#define kc25   274.6   // kinetic coef for CO2 at 25 C, microbars  
#define ko25   419.8   // kinetic coef for O2 at 25C,  millibars 
#define o2     210.0   // oxygen concentration  mmol mol-1  
 

        // tau is computed on the basis of the Specificity factor (102.33) 
        // times Kco2/Kh2o (28.38) to convert for value in solution
        // to that based in air/
        // The old value was 2321.1.

        // New value for Quercus robor from Balaguer et al. 1996
        // Similar number from Dreyer et al. 2001, Tree Physiol, tau= 2710

#define  tau25  2904.12    //  tau coefficient
        //  Arrhenius constants
        //  Eact for Michaelis-Menten const. for KC, KO and dark respiration
        //  These values are from Harley
#define  ekc   80500.0     // Activation energy for K of CO2; J mol-1  
#define  eko   14500.0     // Activation energy for K of O2, J mol-1
#define  erd   38000.0     // activation energy for dark respiration, eg Q10=2 
#define  ektau -29000.0  // J mol-1 (Jordan and Ogren, 1984)
#define tk_25  298.16    // absolute temperature at 25 C
#define toptvc 301.0    // optimum temperature for maximum carboxylation
#define toptjm 301.0    // optimum temperature for maximum electron transport
#define eabole 45162      // activation energy for bole respiration for Q10 = 2.02


        // Constants for leaf energy balance

#define  sigma   5.67e-08   // Stefan-Boltzmann constant W M-2 K-4 
#define  cp      1005.         // Specific heat of air, J KG-1 K-1 
#define mass_air 29.0     // Molecular weight of air, g mole-1 
#define mass_CO2 44.0       // molecular weight of CO2, g mole-1
#define dldt     -2370.0      // Derivative of the latent heat of vaporization
        
#define ep        0.98                     // emissivity of leaves
#define epm1      0.02                    // 1- ep  
#define epsoil    0.98               // Emissivity of soil  
#define epsigma   5.5566e-8           // ep*sigma  
#define epsigma2  11.1132e-8       // 2*ep*sigma
#define epsigma4  22.2264e-8       //  4.0 * ep * sigma
#define epsigma6  33.3396e-8       //  6.0 * ep * sigma 
#define epsigma8  44.448e-8        //  8.0 * ep * sigma 
#define epsigma12 66.6792e-8       // 12.0 * ep * sigma
#define betfact   1.5                 // multiplication factor for aerodynamic 
                                                  // sheltering, based on work by Grace and Wilson 
        //  constants for the polynomial equation for saturation vapor pressure-T function, es=f(t)
#define a1en      617.4
#define a2en      42.22
#define a3en      1.675
#define a4en      0.01408
#define a5en      0.0005818

             
  
        // Minimum stomatal resistance, s m-1.  
#define rsm       145.0
#define brs       60.0      // curvature coeffient for light response

        //   leaf quantum yield, electrons  
#define  qalpha   0.22
#define  qalpha2  0.0484   // qalpha squared, qalpha2 = pow(qalpha, 2.0);

        // Leaf dimension. geometric mean of length and width (m)
#define  lleaf    0.1       // leaf length, m

      
        // Diffusivity values for 273 K and 1013 mb (STP) using values from Massman (1998) Atmos Environment
        // These values are for diffusion in air.  When used these values must be adjusted for
        // temperature and pressure
        // nu, Molecular viscosity 
                

#define  nuvisc  13.27        // mm2 s-1
#define  nnu     0.00001327  // m2 s-1

        // Diffusivity of CO2 

#define  dc    13.81         // mm2 s-1
#define  ddc   0.00001381    // m2 s-1

        //   Diffusivity of heat

#define  dh    18.69;         // mm2 s-1
#define  ddh   0.00001869    // m2 s-1
            

        //  Diffusivity of water vapor 

#define  dv   21.78        // mm2 s-1
#define  ddv   0.00002178    // m2 s-1 
 




