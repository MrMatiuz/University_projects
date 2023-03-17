/*  2D axe-simm  */

/*
    CONTENTS
    void main( void );
    void inp_param(void);
    void memo_on(void);
    void inp_reg(void);
    double step(void);
    void new_gdf(int k,int l);

    void deriv0(void);
    void deriv1(void);
    void deriv(int k);
    void t_deriv(int k, int l);
    double cnvt(int j, int k, int l);
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "func2d.h"

/* #include <time.h> */

 double *gdf[5];
 double AX0,AX1,AY0,AY1 ;
 double dx,dy,dt;
 int kt,lt;
 int MAX_STP;
 double MAX_T,t0,t1,HZ,COORD;
 double fly0[5],fly1[5],*flx0[5],*flx1[5];
 double *d_dx0[5],*d_dx1[5],*d_dy0[5],*d_dy1[5],*d_dt0[5],*d_dt1[5];
 double gdi[5];
 double rqp[3][2],r_o,u_o,p_o,om,s1,s2,s3,g[2];
 FILE *fpr,*fpw;
 double GAM;
 
/*--------------------------------------------------------------------*/
     void main( void )
     {
     int nstp=0;
     float dt1;
     inp_param();
     memo_on();
     inp_reg();
     while( (nstp < MAX_STP)&&(t1 < MAX_T) )
         {
         t1=t0+dt;
         dt1=step();
         dt=HZ*dt1;
         t0=t1;  nstp++;
         printf("STP=%d  T=%e  dt=%e \n",nstp,t0,dt);
         }
     out_reg();
     }
/*--------------------------------------------------------------------*/
    void inp_param(void)
    {
    int i;
    fpr=fopen("par_i.dat","r");
    fscanf(fpr,"%d %d",&kt,&lt);
    fscanf(fpr,"%d %lg %lg %lg",&MAX_STP,&MAX_T,&HZ,&COORD);
    fscanf(fpr,"%lg %lg",&AX0,&AX1);
    fscanf(fpr,"%lg %lg",&AY0,&AY1);
    printf("From inp_param: kt=%d  lt=%d \n",kt,lt);
    for(i=0; i < 5; i++)
       fscanf(fpr,"%lg",&gdi[i]);
    dy=(AY1-AY0)/(float)lt;
    dx=(AX1-AX0)/(float)kt;
    fclose(fpr);
    }
/*--------------------------------------------------------------------*/
    void memo_on(void)
    {
    int i;
    for(i=0; i < 5; i++)
      {
      gdf[i]=(double*)malloc(kt*lt*sizeof(double));
      }

    for(i=0; i < 5; i++)
      {
      flx0[i]=(double*)malloc((lt+1)*sizeof(double));
      flx1[i]=(double*)malloc((lt+1)*sizeof(double));
      d_dx0[i]=(double*)malloc((lt+1)*sizeof(double));
      d_dx1[i]=(double*)malloc((lt+1)*sizeof(double));
      d_dy0[i]=(double*)malloc((lt+1)*sizeof(double));
      d_dy1[i]=(double*)malloc((lt+1)*sizeof(double));
      d_dt0[i]=(double*)malloc((lt+1)*sizeof(double));
      d_dt1[i]=(double*)malloc((lt+1)*sizeof(double));
      }
      
    }
/*--------------------------------------------------------------------*/
    void inp_reg(void)
    {
    int i,kt0,lt0,len_gdf;
    
    fpr=fopen("reg_i.dat","rb");
    fread(&kt0,sizeof(int),1,fpr);
    fread(&lt0,sizeof(int),1,fpr);
    printf("From inp_reg:  kt=%d   lt=%d \n",kt,lt);
    printf("From inp_reg: kt0=%d  lt0=%d \n",kt0,lt0);
    if( (kt!=kt0)||(lt!=lt0) )
      {
      printf("ATT: wrong grid-param: \n");
      printf("     kt=%d   kt0=%d \n",kt,kt0);
      printf("     lt=%d   lt0=%d \n",lt,lt0);
      exit(-1);
      }
    len_gdf=kt*lt;
    fread(&t0,sizeof(double),1,fpr);
    fread(&dt,sizeof(double),1,fpr);
    for(i=0; i < 5; i++)
      fread(gdf[i],sizeof(double),len_gdf,fpr);
    fclose(fpr);
    }
/*--------------------------------------------------------------------*/
    double step(void)
    {
    int i,k,l;
    double dtx,dty,tmp,dt0;
    
    dtx=777777.;  dty=777777.;  dt0=777777.;
    
/*  x-flows for k=0  */
    deriv0();
    tmp=x_flow0();  dtx=min(tmp,dtx);
       
/*  loop: k=0,...kt-1   */

    for(k=0; k < kt; k++)
      {
      
/* x,y,t- derivatives */
      if( k < kt-2 )
                    deriv(k+1);
      else          deriv1();

/*  x-flows  */
      if( k < kt-1 ) { tmp=x_flowk(k); dtx=min(tmp,dtx); }
      else           x_flow1();
      
/*  y-flow for l=0  */
      tmp=y_flow0(k);  dty=min(tmp,dty);
/*  loop:  l=0,1,...lt-1  */
                          
      for(l=0; l < lt; l++)
        {
/*  y-flow for this l  */
        if( l < lt-1)
             {  tmp=y_flow(k,l);  dty=min(tmp,dty); }
        else {  tmp=y_flow1(k);   dty=min(tmp,dty); }
/*  new gdf and step-size  */
        tmp=dty+dtx;
        tmp=dty*dtx/tmp;
        dt0=min(tmp,dt0);
        new_gdf(k,l);
        memcpy(fly0,fly1,5*sizeof(double));
        }
/* end l_loop  */
        for(i=0; i < 5; i++)
          {
          memcpy(flx0[i],flx1[i],lt*sizeof(double));
          memcpy(d_dx0[i],d_dx1[i],lt*sizeof(double));
          memcpy(d_dy0[i],d_dy1[i],lt*sizeof(double));
          memcpy(d_dt0[i],d_dt1[i],lt*sizeof(double));
          }
      }
/*  end k_loop  */
    return dt0;
    }
/*--------------------------------------------------------------------*/
    void new_gdf(int k,int l)
    {
    static double rhs[5],alf[5];
    double yy05,rx,ry,tmp,qr,e;
    int i,j;
    
    i=IGD(k,l);
    GAM=G_A_M(gdf[4][i]);
    yy05=AY0+(l+0.5)*dy;
    rx=dt/dx;  ry=dt/dy;
    e=gdf[3][i]/(GAM-1.)+0.5*gdf[0][i]*(gdf[1][i]*gdf[1][i]+
                                        gdf[2][i]*gdf[2][i]);
    rhs[0]=gdf[0][i]*gdf[2][i];
    rhs[1]=rhs[0]*gdf[1][i];
    rhs[2]=rhs[0]*gdf[2][i];
    rhs[3]=(e+gdf[3][i])*gdf[2][i];
//ATT Only plosk for concentration
    rhs[4]=0.; 
    alf[0]=gdf[0][i];
    alf[1]=gdf[0][i]*gdf[1][i];
    alf[2]=gdf[0][i]*gdf[2][i];
    alf[3]=e;
    alf[4]=gdf[0][i]*gdf[4][i];
    for(j=0; j < 5; j++)
      {
      tmp=rx*(flx1[j][l]-flx0[j][l])+ry*(fly1[j]-fly0[j]);
      alf[j]+=-tmp-dt*COORD*rhs[j]/yy05;
      }
    gdf[0][i]=alf[0];
    tmp=1./gdf[0][i];
    gdf[1][i]=alf[1]*tmp;  gdf[2][i]=alf[2]*tmp;
    qr=0.5*(alf[1]*alf[1]+alf[2]*alf[2])*tmp;
    gdf[3][i]=(GAM-1.)*(alf[3]-qr);
    gdf[4][i]=alf[4]*tmp;
    }
/*--------------------------------------------------------------------*/
    void deriv0(void)
    {
    int i,l;
    for(l=0; l < lt; l++)
      for(i=0; i < 5; i++)
       {
       d_dx0[i][l]=0.;
       d_dy0[i][l]=0.;
       d_dt0[i][l]=0.;
       }
    }
/*--------------------------------------------------------------------*/
    void deriv1(void)
    {
    int i,l;
    for(l=0; l < lt; l++)
      for(i=0; i < 5; i++)
       {
       d_dx1[i][l]=0.;
       d_dy1[i][l]=0.;
       d_dt1[i][l]=0.;
       }
    }
/*--------------------------------------------------------------------*/
    void deriv(int k)
    {
    double tmp1,tmp2;
    int i,l,n,n1,n2;

    for(l=0; l < lt; l++)
     {
/* x,y deriv */
     for(i=0; i < 5; i++)
      {
      n=IGD(k,l);
      if(l==0 || l==lt-1) 
                          d_dy1[i][l]=0.;
      else
          {
          n1=IGD(k,l-1);  n2=IGD(k,l+1);
          tmp1=(gdf[i][n]-gdf[i][n1])/dy;
          tmp2=(gdf[i][n2]-gdf[i][n])/dy;
          d_dy1[i][l]=lim(tmp1,tmp2);
          }
      n1=IGD(k-1,l);  n2=IGD(k+1,l);
      tmp1=(gdf[i][n]-gdf[i][n1])/dx;
      tmp2=(gdf[i][n2]-gdf[i][n])/dx;
      d_dx1[i][l]=lim(tmp1,tmp2);
/* end i-loop */
      }
/* t_deriv */
      t_deriv(k,l);     
/* end l-loop */
     }

    } 
/*--------------------------------------------------------------------*/
   void t_deriv(int k, int l)
   {
   int i;
   double div,tmp,c2;

   div=d_dx1[1][l]+d_dy1[2][l];
   i=IGD(k,l);
   GAM=G_A_M(gdf[4][i]);
   c2=GAM*gdf[3][i]/gdf[0][i];
   d_dt1[0][l]=-cnvt(0,k,l)-gdf[0][i]*div;
   d_dt1[1][l]=-cnvt(1,k,l)-d_dx1[3][l]/gdf[0][i];
   d_dt1[2][l]=-cnvt(2,k,l)-d_dy1[3][l]/gdf[0][i];
   d_dt1[3][l]=-cnvt(3,k,l)-c2*gdf[0][i]*div;
   d_dt1[4][l]=-cnvt(4,k,l);
   }
/*--------------------------------------------------------------------*/
   double cnvt(int j, int k, int l)
   {
   int i;
   i=IGD(k,l);   
   return (gdf[1][i]*d_dx1[j][l]+gdf[2][i]*d_dy1[j][l]);
   }
/*--------------------------------------------------------------------*/
