/*************************************************************************
  readcontrolf.c
  --------
  To Read a model control file.  
  A control file contains input/output files' path and names.
*****************************************************************************/
#include "beps.h"

void readconf(conf)
char conf[];		/* control file name */
{

	char field_descr[255];	/* Field description in control file */
	FILE *fcon;		/* Control file pointer */

/*      If can not open control parameter file */
        if ((fcon = fopen(conf,"r")) == NULL) 
		{
          printf("Can not open parameter file <%s>, exitting program ...\n",conf);
		  exit(0);
        }

/*	Read control file parameters */
    fscanf(fcon,"%s %d", field_descr, &yr);
    fscanf(fcon,"%s %d", field_descr, &bgn_day);
    fscanf(fcon,"%s %d", field_descr, &end_day);
	fscanf(fcon,"%s %d", field_descr, &nlines);
	fscanf(fcon,"%s %d", field_descr, &npixels);
	fscanf(fcon,"%s %s", field_descr, lc_fn);
	fscanf(fcon,"%s %s", field_descr, lai_fn);
	fscanf(fcon,"%s %s", field_descr, stxt_fn);
	fscanf(fcon,"%s %s", field_descr, ci_fn);
	fscanf(fcon,"%s %s", field_descr, lon_fn);
	fscanf(fcon,"%s %s", field_descr, lat_fn);
	fscanf(fcon,"%s %s", field_descr, st_fn);
	fscanf(fcon,"%s %s", field_descr, sw_fn);
	fscanf(fcon,"%s %s", field_descr, sdp_fn);
	fscanf(fcon,"%s %s", field_descr, r_fn);
	fscanf(fcon,"%s %s", field_descr, t_fn);
	fscanf(fcon,"%s %s", field_descr, h_fn);
	fscanf(fcon,"%s %s", field_descr, p_fn);
	fscanf(fcon,"%s %s", field_descr, wd_fn);

	fscanf(fcon,"%s %s", field_descr, laiy_fn);
	fscanf(fcon,"%s %s", field_descr, nppy_fn);
	fscanf(fcon,"%s %s", field_descr, ccd_fn);
	fscanf(fcon,"%s %s", field_descr, cssd_fn);
	fscanf(fcon,"%s %s", field_descr, csmd_fn);
	fscanf(fcon,"%s %s", field_descr, cfsd_fn);
	fscanf(fcon,"%s %s", field_descr, cfmd_fn);
	fscanf(fcon,"%s %s", field_descr, csm_fn);
	fscanf(fcon,"%s %s", field_descr, cm_fn);
	fscanf(fcon,"%s %s", field_descr, cs_fn);
	fscanf(fcon,"%s %s", field_descr, cp_fn);

	fscanf(fcon,"%s %s", field_descr, fp4outp1);
	fscanf(fcon,"%s %s", field_descr, fp4outp2);
	fscanf(fcon,"%s %s", field_descr, fp4outp3);

	
/*	Output control parameters */
	printf("\n Simulating...\n");
	printf(" day: %d - %d   year: %d\n",bgn_day,end_day,yr);
//	printf(" LC: %s CI: %s  soilt: %s\n",lc_fn,ci_fn,stxt_fn);
//	printf(" output: %s  %s  %s\n",fp4outp1,fp4outp2,fp4outp3);


fclose(fcon);

}



