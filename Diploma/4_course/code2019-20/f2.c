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


extern double *gdf[6],*x,*y,*fllr[5],*fldu[5],*xc,*yc,
       *nlr[3],*tlr[2],*ndu[3],*tdu[2],*vol;
extern double *dxflr[5],*dyflr[5],*dxfdu[5],*dyfdu[5];
extern double t,dt,t_max,coord,turb,Mach,Re_ref,u_ref,mu_ref;
extern double lxmin,lymin,mu_max,H_Z;
extern int n_stp,max_stp,inp_par,IT,JT;

extern double *wlr[2],*wdu[2],*gdf_[6];
extern double tau;
extern int i_qgd;

extern double *mu_turb_lr, *mu_turb_du, *mu_To_du, *mu_Ti_du, *mu_To_lr, *mu_Ti_lr;


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
 fprintf(fpw," VARIABLES = \"X\", \"Y\", \"D\", \"U\", \"V\", \"P\",\"MU_TURB\",\"mu_turb_lr\",\"mu_turb_du\", \"Ti_du\", \"To_du\", \"Ti_lr\", \"To_lr\"   \n");
 fprintf(fpw," ZONE F=BLOCK  I=%d  J=%d \n", JT-3, IT-3);
 fprintf(fpw," VARLOCATION=([3-13]=CELLCENTERED) \n");

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
    fprintf(fpw," %12.4e ", 0.25 * (mu_turb_lr[Ilr(i,j)] + mu_turb_lr[Ilr(i+1,j)] + mu_turb_du[Idu(i,j)] + mu_turb_du[Idu(i,j+1)])); // mu_turb_arr[IC(i,j)]
    iprn++;
    if(iprn > 6) { iprn=1; fprintf(fpw,"\n"); }
  }
 fprintf(fpw,"\n");

 iprn=1;
 for(i=2; i < IT-2; i++)
  for(j=2; j < JT-2; j++)
   {
    fprintf(fpw," %12.4e ", mu_turb_lr[Ilr(i,j)]);
    iprn++;
    if(iprn > 6) { iprn=1; fprintf(fpw,"\n"); }
   }
 fprintf(fpw,"\n");

 iprn=1;
 for(i=2; i < IT-2; i++)
  for(j=2; j < JT-2; j++)
   {
    fprintf(fpw," %12.4e ", mu_turb_du[Idu(i,j)]);
    iprn++;
    if(iprn > 6) { iprn=1; fprintf(fpw,"\n"); }
   }
 fprintf(fpw,"\n");

 iprn=1;
 for(i=2; i < IT-2; i++)
  for(j=2; j < JT-2; j++)
   {
    fprintf(fpw," %12.4e ", mu_Ti_du[Idu(i,j)]);
    iprn++;
    if(iprn > 6) { iprn=1; fprintf(fpw,"\n"); }
   }
 fprintf(fpw,"\n");

 iprn=1;
 for(i=2; i < IT-2; i++)
  for(j=2; j < JT-2; j++)
   {
    fprintf(fpw," %12.4e ", mu_To_du[Idu(i,j)]);
    iprn++;
    if(iprn > 6) { iprn=1; fprintf(fpw,"\n"); }
   }
 fprintf(fpw,"\n");

  iprn=1;
 for(i=2; i < IT-2; i++)
  for(j=2; j < JT-2; j++)
   {
    fprintf(fpw," %12.4e ", mu_Ti_lr[Ilr(i,j)]);
    iprn++;
    if(iprn > 6) { iprn=1; fprintf(fpw,"\n"); }
   }
 fprintf(fpw,"\n");

 iprn=1;
 for(i=2; i < IT-2; i++)
  for(j=2; j < JT-2; j++)
   {
    fprintf(fpw," %12.4e ", mu_To_lr[Ilr(i,j)]);
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

 mu_turb_lr=(double*)malloc(len*sizeof(double));
 mu_turb_du=(double*)malloc(len*sizeof(double));
 mu_To_du=(double*)malloc(len*sizeof(double));
 mu_Ti_du=(double*)malloc(len*sizeof(double));
 mu_To_lr=(double*)malloc(len*sizeof(double));
 mu_Ti_lr=(double*)malloc(len*sizeof(double));

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
  dxflr[k]=(double*)malloc(len*sizeof(double));
  dyflr[k]=(double*)malloc(len*sizeof(double));
  fldu[k]=(double*)malloc(len*sizeof(double));
  dxfdu[k]=(double*)malloc(len*sizeof(double));
  dyfdu[k]=(double*)malloc(len*sizeof(double));
  } 
 
 for(k=0; k < 6; k++)
  {
  gdf[k]=(double*)malloc(len*sizeof(double));
  gdf_[k]=(double*)malloc(len*sizeof(double));
  }
for(k=0; k < 2; k++)
  {
  wlr[k]=(double*)malloc(len*sizeof(double));
  wdu[k]=(double*)malloc(len*sizeof(double));
  }
}
