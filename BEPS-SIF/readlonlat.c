/*************************************************************************
    Read longitude and latitude from a file in a given line.
*****************************************************************************/
#include <malloc.h>
#include <stdio.h>
#include <string.h>
#include "beps.h"

void readlonlat(lon,lat)
float **lon,**lat; 
{
	short line,pix;
	float lcv;
	FILE *lon_ptr,*lat_ptr;
	
	/*      Open longitude file */
        if ((lon_ptr=fopen(lon_fn, "r")) == NULL) 
	{ printf("\n Unable to open lon. file, exitting...\n");
          exit(0);
        }
        /*      Open latitude file */
        if ((lat_ptr=fopen(lat_fn, "r")) == NULL) 
	{ printf("\n Unable to open lat. file, exitting...\n");
          exit(0);
        }
	
	for (line = 0; line <nlines; line++) 
		for(pix=0; pix<npixels; pix++) 
		{  
                 fscanf(lon_ptr,"%f ",&lcv);
                 lon[line][pix]=lcv;
                 fscanf(lat_ptr,"%f ",&lcv);
                 lat[line][pix]=lcv;

		}
	
	return;
}
