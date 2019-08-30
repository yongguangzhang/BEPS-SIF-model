/*************************************************************************
  Description:
  -----------
  To read avg. LAI, ann_NPP and carbon pools data for a pixel
*****************************************************************************/

#include"beps.h"

void 
readinput2(lai_y,nppy,Ccd,Cssd,Csmd,Cfsd,Cfmd,Csm,Cm,Cs,Cp)
float **lai_y,**nppy; 
float **Ccd,**Cssd,**Csmd,**Cfsd,**Cfmd,**Csm,**Cm,**Cs,**Cp;
{
	short line,pix;
	float laiv,nppv_y;
	float ccdv,cssdv,csmdv,cfsdv,cfmdv,csmv,cmv,csv,cpv; 
	
	FILE *laiy_ptr,*nppy_ptr;
	FILE *ccdf_ptr,*cssdf_ptr,*csmdf_ptr,*cfsdf_ptr,*cfmdf_ptr,*csmf_ptr,*cmf_ptr,*csf_ptr,*cpf_ptr;		

	/*      Open lai yr avg. file */
  	if ((laiy_ptr=fopen(laiy_fn, "r")) == NULL) {
          printf("\n Unable to open LC file, exitting program...\n\n");
	  exit(0);
    	}
	/*      Open NPP year avg. file */
	if ((nppy_ptr=fopen(nppy_fn, "rb")) == NULL) {
          printf("\n Unable to open yrnpp file,  exitting program...\n\n");
	  exit(0);
    	}

	  for (line = 0; line <nlines; line++) 
		for(pix=0; pix<npixels; pix++) 
		{  
		 fscanf(laiy_ptr,"%f ",&laiv);
		 lai_y[line][pix]=laiv;
		 fread(&nppv_y,sizeof(float),1,nppy_ptr);
		 nppy[line][pix]=nppv_y;
		}
	fclose(laiy_ptr);
	fclose(nppy_ptr);
                

    /*      Open soil carbon  and soil int. data files */
  	if ((ccdf_ptr=fopen(ccd_fn, "r")) == NULL) {
          printf("\n Unable to open file 5,  exitting program ...\n\n");
	  exit(0);
    	}
  	if ((cssdf_ptr=fopen(cssd_fn, "r")) == NULL) {
          printf("\n Unable to open file 6,  exitting program ...\n\n");
	  exit(0);
    	}
  	if ((csmdf_ptr=fopen(csmd_fn, "r")) == NULL) {
          printf("\n Unable to open file 7,  exitting program ...\n\n");
	  exit(0);
    	}
  	if ((cfsdf_ptr=fopen(cfsd_fn, "r")) == NULL) {
          printf("\n Unable to open file 8,  exitting program ...\n\n");
	  exit(0);
    	}
  	if ((cfmdf_ptr=fopen(cfmd_fn, "r")) == NULL) {
          printf("\n Unable to open file 9,  exitting program ...\n\n");
	  exit(0);
    	}
  	if ((csmf_ptr=fopen(csm_fn, "r")) == NULL) {
          printf("\n Unable to open file 10,  exitting program ...\n\n");
	  exit(0);
    	}
  	if ((cmf_ptr=fopen(cm_fn, "r")) == NULL) {
          printf("\n Unable to open file 11,  exitting program ...\n\n");
	  exit(0);
    	}
  	if ((csf_ptr=fopen(cs_fn, "r")) == NULL) {
          printf("\n Unable to open file 12,  exitting program ...\n\n");
	  exit(0);
    	}
  	if ((cpf_ptr=fopen(cp_fn, "r")) == NULL) {
          printf("\n Unable to open file 13,  exitting program ...\n\n");
	  exit(0);
    	}

/*	  Read carbon data for each pix */
	  for (line = 0; line <nlines; line++) 
		for(pix=0; pix<npixels; pix++) 
		{  
			fscanf(ccdf_ptr,"%f ",&ccdv);
			Ccd[line][pix]=ccdv;
			fscanf(cssdf_ptr,"%f ",&cssdv);
			Cssd[line][pix]=cssdv;
			fscanf(csmdf_ptr,"%f ",&csmdv);
			Csmd[line][pix]=csmdv;
			fscanf(cfsdf_ptr,"%f ",&cfsdv);
			Cfsd[line][pix]=cfsdv;
			fscanf(cfmdf_ptr,"%f ",&cfmdv);
			Cfmd[line][pix]=cfmdv;
			fscanf(csmf_ptr,"%f ",&csmv);
			Csm[line][pix]=csmv;
			fscanf(cmf_ptr,"%f ",&cmv);
			Cm[line][pix]=cmv;
			fscanf(csf_ptr,"%f ",&csv);
			Cs[line][pix]=csv;
			fscanf(cpf_ptr,"%f ",&cpv);
			Cp[line][pix]=cpv;
		}

	fclose(ccdf_ptr);
	fclose(cssdf_ptr);
	fclose(csmdf_ptr);
	fclose(cfsdf_ptr);
	fclose(cfmdf_ptr);
	fclose(csmf_ptr);
	fclose(cmf_ptr);
	fclose(csf_ptr);
	fclose(cpf_ptr);

}
