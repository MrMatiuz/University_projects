/* CONTENTS
void out(void);
void out_tec(void);
void memo_on(void);
*/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "func.h"

#define IN(i,j)  ( (i)*(JT+1)+(j) )
#define IC(i,j)  ( (i)*JT+(j) )
#define Ilr(i,j)  ( (i)*JT+(j) )
#define Idu(i,j)  ( (i)*(JT+1)+(j) )

#define  GAM  1.4


extern double *gdf[6],*gdf_[6],*w[2],*x,*y,*xc,*yc,*fllr[5],*fldu[5],
       *nlr[3],*tlr[2],*ndu[3],*tdu[2],*vol;
extern double *dxf[5],*dyf[5];
extern double t,dt,t_max,coord,turb,Mach,Re_ref,u_ref,mu_ref;
extern double lxmin,lymin,mu_max,H_Z;
extern int n_stp,max_stp,inp_par,IT,JT;


//*****************************************************
void out(void)
{
 FILE *fpw;
 int i,j,k;
  fpw=fopen("reg_o.dat","wb");
  fwrite(&IT,sizeof(int),1,fpw);
  fwrite(&JT,sizeof(int),1,fpw);
  fwrite(&t,sizeof(double),1,fpw);
  fwrite(&dt,sizeof(double),1,fpw);
   for(i=2; i < IT-2; i++)
    for(j=2; j < JT-2; j++)
     for(k=0; k < 5; k++)
      fwrite(&gdf[k][IC(i,j)],sizeof(double),1,fpw);
   fclose(fpw);
}

//*****************************************************
void out_tec(void)
{
 FILE *fpw;
 int i,j,k,iprn;
 fpw=fopen("fld.dat","w");

 fprintf(fpw," TITLE = \" NS2D-data \"  \n");
 fprintf(fpw," VARIABLES = \"X\", \"Y\", \"D\", \"U\", \"V\", \"P\", \"W0\", \"W1\"  \n");
 fprintf(fpw," ZONE F=BLOCK  I=%d  J=%d \n",JT-3,IT-3);
 fprintf(fpw," VARLOCATION=([3-8]=CELLCENTERED) \n");

 iprn=1;
 for(i=2; i < IT-1; i++)
  for(j=2; j < JT-1; j++)
   {
    fprintf(fpw," %12.4e ",x[IN(i,j)]);
    iprn++;
    if(iprn > 6) { iprn=1; fprintf(fpw,"\n"); }
   }
 fprintf(fpw,"\n");

 iprn=1;
 for(i=2; i < IT-1; i++)
  for(j=2; j < JT-1; j++)
   {
    fprintf(fpw," %12.4e ",y[IN(i,j)]);
    iprn++;
    if(iprn > 6) { iprn=1; fprintf(fpw,"\n"); }
   }
 fprintf(fpw,"\n");

 iprn=1;
 for(i=2; i < IT-2; i++)
  for(j=2; j < JT-2; j++)
   {
    fprintf(fpw," %12.4e ",gdf[0][IC(i,j)]);
    iprn++;
    if(iprn > 6) { iprn=1; fprintf(fpw,"\n"); }
   }
 fprintf(fpw,"\n");

 iprn=1;
 for(i=2; i < IT-2; i++)
  for(j=2; j < JT-2; j++)
   {
    fprintf(fpw," %12.4e ",gdf[1][IC(i,j)]);
    iprn++;
    if(iprn > 6) { iprn=1; fprintf(fpw,"\n"); }
   }
 fprintf(fpw,"\n");

 iprn=1;
 for(i=2; i < IT-2; i++)
  for(j=2; j < JT-2; j++)
   {
    fprintf(fpw," %12.4e ",gdf[2][IC(i,j)]);
    iprn++;
    if(iprn > 6) { iprn=1; fprintf(fpw,"\n"); }
   }
 fprintf(fpw,"\n");

 iprn=1;
 for(i=2; i < IT-2; i++)
  for(j=2; j < JT-2; j++)
   {
    fprintf(fpw," %12.4e ",gdf[3][IC(i,j)]);
    iprn++;
    if(iprn > 6) { iprn=1; fprintf(fpw,"\n"); }
   }
 fprintf(fpw,"\n");
 iprn=1;
 for(i=2; i < IT-2; i++)
  for(j=2; j < JT-2; j++)
   {
    fprintf(fpw," %12.4e ",w[0][IC(i,j)]);
    iprn++;
    if(iprn > 6) { iprn=1; fprintf(fpw,"\n"); }
   }
 fprintf(fpw,"\n");
 iprn=1;
 for(i=2; i < IT-2; i++)
  for(j=2; j < JT-2; j++)
   {
    fprintf(fpw," %12.4e ",w[1][IC(i,j)]);
    iprn++;
    if(iprn > 6) { iprn=1; fprintf(fpw,"\n"); }
   }
 fprintf(fpw,"\n");

 fclose(fpw);
}

//*****************************************************
void memo_on(void)
{
 int k,len;
 len=(IT+1)*(JT+1);

 x=(double*)malloc(len*sizeof(double));
 y=(double*)malloc(len*sizeof(double));
 vol=(double*)malloc(len*sizeof(double));
 xc=(double*)malloc(len*sizeof(double));
 yc=(double*)malloc(len*sizeof(double));

 for(k=0; k < 2; k++)
  {
  tlr[k]=(double*)malloc(len*sizeof(double));
  tdu[k]=(double*)malloc(len*sizeof(double));
  }

 for(k=0; k < 3; k++)
  {
  nlr[k]=(double*)malloc(len*sizeof(double));
  ndu[k]=(double*)malloc(len*sizeof(double));
  }

 for(k=0; k < 5; k++)
  {
  fllr[k]=(double*)malloc(len*sizeof(double));
  dxf[k]=(double*)malloc(len*sizeof(double));
  dyf[k]=(double*)malloc(len*sizeof(double));
  fldu[k]=(double*)malloc(len*sizeof(double));
  }

 for(k=0; k < 6; k++)
  {
  gdf[k]=(double*)malloc(len*sizeof(double));
  gdf_[k]=(double*)malloc(len*sizeof(double));
  }

for(k=0; k < 2; k++)
  w[k]=(double*)malloc(len*sizeof(double));
}
