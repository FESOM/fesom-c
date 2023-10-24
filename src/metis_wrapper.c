
/* SF 06/2001 */
/* NR 07/2014, 12/2015 */
/* The partinioner Interface to METIS-5.1 */

/* INPUT VARIABLES */
/*   n:     dimension of system (number of equations=number of variables) */
/*   ptr :  rowptr or colptr (compressed sparse matrix)  */
/*   adj :  colind or rowind (compressed sparse matrix)  */
/*   np  :  partitioning for np processors */

/* OUTPUT VARIABLES */
/*   part : partitioning vector */
/*   ec :   edge cut  */
/*   ierr:  Metis error code */  

#include "metis.h" 
#include <stdio.h>


void metis_wrapper(idx_t *n, idx_t *ptr, idx_t *adj, idx_t *np, idx_t *part, idx_t *ec, int *ierr)
{

int i ;
idx_t opt[METIS_NOPTIONS];  
idx_t ncon;


if (*np==1) { for(i=0;i<*n;i++) part[i]=0; return;}

 ncon=1;
 METIS_SetDefaultOptions(opt);

 opt[METIS_OPTION_CONTIG]    = 0;  /* 1: contiguous partitions, please */
                                   /* With more weights, this makes the partitions uglier... */
                                   /* Ignored by METIS_PartGraphRecursive */

 opt[METIS_OPTION_NCUTS]   = 10; /* Build NCUTS partitions, choose the best */
 opt[METIS_OPTION_NITER]   = 40;/* higher => better quality, more compute time. Default: 10 */


 /* Try, which one gives better partitioning. Often, recursive is better for load balancing, */
    /* Kway for edge cut */
  *ierr = METIS_PartGraphRecursive(n,&ncon,ptr,adj,NULL,NULL,NULL,np,NULL,NULL,opt,ec,part);
  // *ierr = METIS_PartGraphKway(n,&ncon,ptr,adj,NULL,NULL,NULL,np,NULL,NULL,opt,ec,part);

}
