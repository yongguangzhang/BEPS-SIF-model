/*************************************************************************
  To read constant input data for a pixel
*****************************************************************************/

#include"beps.h"

void 
readinput1(lc,ci,stxt,st2,sw2,sdp2)

float **lc,**ci,**stxt; 
float **st2,**sw2,**sdp2;

{
	short line,pix;
	int lcv,stxtv;
	float stv,swv,sdpv,civ;
	FILE *lcf_ptr,*ci_ptr,*stxt_ptr;
	FILE *st_ptr,*sw_ptr,*sdp_ptr;

	/*      Open landcover file */
  	if ((lcf_ptr=fopen(lc_fn, "r")) == NULL) {
          printf("\n Unable to open lc file, exitting program...\n");
	  exit(0);
    	}
    /*      Open clumping index file */
    if ((ci_ptr=fopen(ci_fn, "r")) == NULL) {
          printf("\n Unable to open CI file,  exitting program...\n\n");
          exit(0);
        }
    /*      Open soil texture file */
  	if ((stxt_ptr=fopen(stxt_fn, "r")) == NULL) {
          printf("\n Unable to open file 14,  exitting program ...\n\n");
	  exit(0);
    	}

	  for (line = 0; line <nlines; line++) 
		for(pix=0; pix<npixels; pix++) 
		{  
	/*  Read landcover data for each pix */
		 fscanf(lcf_ptr,"%d ",&lcv);
		 lc[line][pix]=(float)lcv;
    
	/*  Read clumping index */
         	fscanf(ci_ptr,"%f ",&civ);
         	ci[line][pix]=(float)(0.01*civ);

	/*  Read soiltxt data for each pix */
		 fscanf(stxt_ptr,"%d ",&stxtv);
		 stxt[line][pix]=(float)stxtv;

		}
	fclose(lcf_ptr);
	fclose(ci_ptr);
	fclose(stxt_ptr);
                 
    /*      Open soil int. data files */
  	if ((st_ptr=fopen(st_fn, "r")) == NULL) {
          printf("\n Unable to open file 15,  exitting program ...\n\n");
	  exit(0);
    	}
  	if ((sw_ptr=fopen(sw_fn, "r")) == NULL) {
          printf("\n Unable to open file 16,  exitting program ...\n\n");
	  exit(0);
    	}
  	if ((sdp_ptr=fopen(sdp_fn, "r")) == NULL) {
          printf("\n Unable to open file 17,  exitting program ...\n\n");
	  exit(0);
    	}
/*	  Read soiltemp, water and snowdepth */
	  for (line = 0; line <nlines; line++) 
		for(pix=0; pix<npixels; pix++) 
		{  
			fscanf(st_ptr,"%f ",&stv);
			st2[line][pix]=stv;
			fscanf(sw_ptr,"%f ",&swv);
			sw2[line][pix]=swv;
			fscanf(sdp_ptr,"%f ",&sdpv);
			sdp2[line][pix]=sdpv;
		}
	fclose(st_ptr);
	fclose(sw_ptr);
	fclose(sdp_ptr);


}
