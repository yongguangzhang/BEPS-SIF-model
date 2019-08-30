/*************************************************************
  To read LAI data 
*************************************************************/

#include"beps.h"

void readlai_all(lai0)
float **lai0; 
{
	int l,p,jday,yd,x;
	float laiv;
	char laif[255];
	FILE *laif_ptr;

    
	for (jday=1; jday<=361; jday=jday+8)
	{
	  yd = yr*1000+jday;

	  sprintf(laif,"%s%d.txt",lai_fn,yd);

  	  if ((laif_ptr=fopen(laif, "r")) == NULL) {
          printf("\n Unable to open laif file ...\n\n");
	  exit(0);   	}

	  /* Read lai data for each pix */
    	  for (l=0; l<nlines; l++) 
		for(p=0; p<npixels; p++) 
		{ 
		  x=l*npixels+p; 
		  fscanf(laif_ptr,"%f ",&laiv);
		  lai0[jday][x]=laiv;
	  	  lai0[jday+1][x]=laiv;
	  	  lai0[jday+2][x]=laiv;
	  	  lai0[jday+3][x]=laiv;
	  	  lai0[jday+4][x]=laiv;
	  	  lai0[jday+5][x]=laiv;
	  	  lai0[jday+6][x]=laiv;
	  	  lai0[jday+7][x]=laiv;
		}

	fclose(laif_ptr);
	}  // end of for
}
