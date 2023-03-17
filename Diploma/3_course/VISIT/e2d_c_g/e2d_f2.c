/*  CONTENTS
    double x_flow0(void);
    void x_flow1(void);
    double x_flowk( int k );
    double  y_flow(int k,int l);
    double  y_flow1(int k);
    double  y_flow0(int k);
    void out_reg(void);
    double lim(double a,double b);
    double riem_c(double c_l,double c_r,double spd);

*/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "func2d.h"
#include "data2d.h"

/*--------------------------------------------------------------------*/
    double x_flow0(void)
    {
    int k,l,i2,j,lom;
    double e;
    static double tmp[5];
    double rb,ub,vb,pb;
    double t_x,t_mp;
    double cl,cr,qi,cb;

    lom=0;  om=0.;   t_x=777777.;
    k=0;
      for(l=0; l < lt; l++)
        {
        i2=IGD(k,l);
        rqp[0][0]=gdi[0];  rqp[0][1]=gdf[0][i2];
        rqp[1][0]=gdi[1];  rqp[1][1]=gdf[1][i2];
        rqp[2][0]=gdi[3];  rqp[2][1]=gdf[3][i2];
        g[0]=G_A_M(gdi[4]); g[1]=G_A_M(gdf[4][i2]);
        if( riemi (lom) )
          {
          printf("ALARM_RIEM FROM x_flowk: k=%d l=%d  \n",k,l);
          printf("L: %e %e %e \n",rqp[0][0],rqp[1][0],rqp[2][0]);
          printf("R: %e %e %e \n",rqp[0][1],rqp[1][1],rqp[2][1]);
          exit(-1);
          }
        rb=r_o; ub=u_o; pb=p_o;
        if( s2 > om )
                    vb=gdi[2];
        else
                    vb=gdf[2][i2];
/*  step t_x  */
        if(s3 > 0.) { t_mp=dx/s3;       t_x=min(t_x,t_mp); }
                     
        cl=gdi[4]; cr=gdf[4][i2];
        qi=ub;
        cb=riem_c(cl,cr,qi);
        GAM=G_A_M(cb);
        e=pb/(GAM-1.)+0.5*rb*(ub*ub+vb*vb);
        tmp[0]=rb*ub;
        tmp[1]=pb+rb*ub*ub;
        tmp[2]=rb*ub*vb;
        tmp[3]=(e+pb)*ub;
        tmp[4]=rb*ub*cb; 
        for(j=0; j < 5; j++)
          flx0[j][l]=tmp[j];
        }
    return t_x;
    }
/*--------------------------------------------------------------------*/
    void x_flow1(void)
    {
    int l,i,n;
    double e;
    static double tmp[5];
      for(l=0; l < lt; l++)
        {
        n=IGD(kt-1,l);
        GAM=G_A_M(gdf[4][n]);
        e=gdf[3][n]/(GAM-1.)+0.5*gdf[0][n]*
          (gdf[1][n]*gdf[1][n]+gdf[2][n]*gdf[2][n]);
        tmp[0]=gdf[0][n]*gdf[1][n];
        tmp[1]=gdf[3][n]+tmp[0]*gdf[1][n];
        tmp[2]=tmp[0]*gdf[2][n];
        tmp[3]=(e+gdf[3][n])*gdf[1][n];
        tmp[4]=tmp[0]*gdf[4][n];
        for(i=0; i < 5; i++)
          flx1[i][l]=tmp[i];
        }
    }
/*--------------------------------------------------------------------*/
    double x_flowk( int k )
    {
    int l,i1,i2,j,lom;
    double e;
    static double tmp[5];
    double rb,ub,vb,pb;
    double t_x,t_mp;
    double cl,cr,qi,cb;
    lom=0;  om=0.;   t_x=777777.;
      for(l=0; l < lt; l++)
        {
        i1=IGD(k,l);    i2=IGD(k+1,l);

        rqp[0][0]=gdf[0][i1]+0.5*dx*d_dx0[0][l]+0.5*dt*d_dt0[0][l]; 
        rqp[1][0]=gdf[1][i1]+0.5*dx*d_dx0[1][l]+0.5*dt*d_dt0[1][l]; 
        rqp[2][0]=gdf[3][i1]+0.5*dx*d_dx0[3][l]+0.5*dt*d_dt0[3][l]; 

        rqp[0][1]=gdf[0][i2]-0.5*dx*d_dx1[0][l]+0.5*dt*d_dt1[0][l]; 
        rqp[1][1]=gdf[1][i2]-0.5*dx*d_dx1[1][l]+0.5*dt*d_dt1[1][l]; 
        rqp[2][1]=gdf[3][i2]-0.5*dx*d_dx1[3][l]+0.5*dt*d_dt1[3][l]; 
        g[0]=G_A_M(gdf[4][i1]); g[1]=G_A_M(gdf[4][i2]);

        if( riemi (lom) )
          {
          printf("ALARM_RIEM FROM x_flowk: k=%d l=%d  \n",k,l);
          printf("L: %e %e %e \n",rqp[0][0],rqp[1][0],rqp[2][0]);
          printf("R: %e %e %e \n",rqp[0][1],rqp[1][1],rqp[2][1]);
          exit(-1);
          }
        rb=r_o; ub=u_o; pb=p_o;
        if( s2 > om )
                    vb=gdf[2][i1];
        else
                    vb=gdf[2][i2];
/*  step t_x  */
        if(s1 < 0.) { t_mp=dx/fabs(s1); t_x=min(t_x,t_mp); }
        if(s3 > 0.) { t_mp=dx/s3;       t_x=min(t_x,t_mp); }
                     
        cl=gdf[4][i1]+0.5*dx*d_dx0[4][l]+0.5*dt*d_dt0[4][l]; 
        cr=gdf[4][i2]-0.5*dx*d_dx1[4][l]+0.5*dt*d_dt1[4][l]; 
        qi=ub;
        cb=riem_c(cl,cr,qi);
        GAM=G_A_M(cb);

        e=pb/(GAM-1.)+0.5*rb*(ub*ub+vb*vb);
        tmp[0]=rb*ub;
        tmp[1]=pb+rb*ub*ub;
        tmp[2]=rb*ub*vb;
        tmp[3]=(e+pb)*ub;
        tmp[4]=rb*ub*cb; 
        for(j=0; j < 5; j++)
          flx1[j][l]=tmp[j];
        }
    return t_x;
    }
/*--------------------------------------------------------------------*/
    double  y_flow(int k,int l)
    {
    int n,i,i1,i2,j,lom;
    double e;
    static double tmp[5];
    double rb,ub,vb,pb;
    double cl,cr,qi,cb;
    double t_y,t_mp;
    lom=0;  t_y=777777.;
    i1=IGD(k,l);  i2=IGD(k,(l+1));

    rqp[0][0]=gdf[0][i1]+0.5*dy*d_dy0[0][l]+0.5*dt*d_dt0[0][l];
    rqp[1][0]=gdf[2][i1]+0.5*dy*d_dy0[2][l]+0.5*dt*d_dt0[2][l];
    rqp[2][0]=gdf[3][i1]+0.5*dy*d_dy0[3][l]+0.5*dt*d_dt0[3][l];

    rqp[0][1]=gdf[0][i2]-0.5*dy*d_dy0[0][l+1]+0.5*dt*d_dt0[0][l+1];
    rqp[1][1]=gdf[2][i2]-0.5*dy*d_dy0[2][l+1]+0.5*dt*d_dt0[2][l+1];
    rqp[2][1]=gdf[3][i2]-0.5*dy*d_dy0[3][l+1]+0.5*dt*d_dt0[3][l+1];
    g[0]=G_A_M(gdf[4][i1]); g[1]=G_A_M(gdf[4][i2]);

    if( riemi (lom) )
      {
      printf("ALARM_RIEM FROM y_flow: k=%d l=%d \n",k,l);
      printf("L: %e %e %e \n",rqp[0][0],rqp[1][0],rqp[2][0]);
      printf("R: %e %e %e \n",rqp[0][1],rqp[1][1],rqp[2][1]);
      exit(-1);
      }
    rb=r_o; vb=u_o; pb=p_o;
    
/*  step t_r  */
    if(s1 < 0.) { t_mp=dy/fabs(s1); t_y=min(t_y,t_mp); }
    if(s3 > 0.) { t_mp=dy/s3;       t_y=min(t_y,t_mp); }
    
    if( s2 > om )
                 ub=gdf[1][i1];
    else
                 ub=gdf[1][i2];

    cl=gdf[4][i1]+0.5*dy*d_dy0[4][l]+0.5*dt*d_dt0[4][l];
    cr=gdf[4][i2]-0.5*dy*d_dy0[4][l+1]+0.5*dt*d_dt0[4][l+1];
    qi=vb;
    cb=riem_c(cl,cr,qi);
    GAM=G_A_M(cb);
    e=pb/(GAM-1.)+0.5*rb*(ub*ub+vb*vb);
    tmp[0]=rb*vb;
    tmp[1]=rb*ub*vb;
    tmp[2]=pb+rb*vb*vb;
    tmp[3]=(e+pb)*vb;
    tmp[4]=cb*rb*vb;
    for(j=0; j < 5; j++)
      fly1[j]=tmp[j];
    return t_y;
    }
/*--------------------------------------------------------------------*/
    double  y_flow1(int k)
    {
    int l,n,i,i1,j,lom;
    double e;
    static double tmp[5];
    double rb,ub,vb,pb;
    double t_y,t_mp;
    l=lt-1;
    lom=0;  t_y=777777.;
    i1=IGD(k,l);
    rqp[0][0]=gdf[0][i1];  rqp[0][1]=gdf[0][i1];
    rqp[1][0]=gdf[2][i1];  rqp[1][1]=-gdf[2][i1];
    rqp[2][0]=gdf[3][i1];  rqp[2][1]=gdf[3][i1];
    g[0]=G_A_M(gdf[4][i1]); g[1]=G_A_M(gdf[4][i1]);

    if( riemi (lom) )
      {
      printf("ALARM_RIEM FROM y_flow1: k=%d l=%d \n",k,l);
      printf("L: %e %e %e \n",rqp[0][0],rqp[1][0],rqp[2][0]);
      printf("R: %e %e %e \n",rqp[0][1],rqp[1][1],rqp[2][1]);
      exit(-1);
      }
    rb=r_o; vb=u_o; pb=p_o;
    
/*  step t_r  */
    if(s1 < 0.) { t_mp=dy/fabs(s1); t_y=min(t_y,t_mp); }
    
    ub=gdf[1][i1];
    GAM=0.5*(g[0]+g[1]);
    e=pb/(GAM-1.)+0.5*rb*(ub*ub+vb*vb);
    tmp[0]=rb*vb;
    tmp[1]=rb*ub*vb;
    tmp[2]=pb+rb*vb*vb;
    tmp[3]=(e+pb)*vb;
    tmp[4]=0.;
    for(j=0; j < 5; j++)
      fly1[j]=tmp[j];
    return t_y;
    }
/*--------------------------------------------------------------------*/
    double  y_flow0(int k)
    {
    int l,n,i,i1,j,lom;
    double e;
    static double tmp[5];
    double rb,ub,vb,pb;
    double t_y,t_mp;
    l=0;
    lom=0;  t_y=777777.;
    i1=IGD(k,l);
    rqp[0][0]= gdf[0][i1];  rqp[0][1]=gdf[0][i1];
    rqp[1][0]=-gdf[2][i1];  rqp[1][1]=gdf[2][i1];
    rqp[2][0]= gdf[3][i1];  rqp[2][1]=gdf[3][i1];
    g[0]=G_A_M(gdf[4][i1]); g[1]=G_A_M(gdf[4][i1]);
    if( riemi (lom) )
      {
      printf("ALARM_RIEM FROM y_flow1: k=%d l=%d \n",k,l);
      printf("L: %e %e %e \n",rqp[0][0],rqp[1][0],rqp[2][0]);
      printf("R: %e %e %e \n",rqp[0][1],rqp[1][1],rqp[2][1]);
      exit(-1);
      }
    rb=r_o; vb=u_o; pb=p_o;
    
/*  step t_r  */
    if(s3 > 0.) { t_mp=dy/s3;       t_y=min(t_y,t_mp); }
    
    ub=gdf[1][i1];
    GAM=0.5*(g[0]+g[1]);
    e=pb/(GAM-1.)+0.5*rb*(ub*ub+vb*vb);
    tmp[0]=rb*vb;
    tmp[1]=rb*ub*vb;
    tmp[2]=pb+rb*vb*vb;
    tmp[3]=(e+pb)*vb;
    tmp[4]=0.;
    for(j=0; j < 5; j++)
      fly0[j]=tmp[j];
    return t_y;
    }
/*--------------------------------------------------------------------*/
    void out_reg(void)
    {
    int i,len_gdf;
    fpw=fopen("reg_o.dat","wb");
    fwrite(&kt,sizeof(int),1,fpw);
    fwrite(&lt,sizeof(int),1,fpw);
    len_gdf=kt*lt;
    fwrite(&t0,sizeof(double),1,fpw);
    fwrite(&dt,sizeof(double),1,fpw);
    for(i=0; i < 5; i++)
      fwrite(gdf[i],sizeof(double),len_gdf,fpw);
    fclose(fpw);
    }
/*--------------------------------------------------------------------*/
double lim(double a,double b)
{
double tmp;
if( a*b <= 0.) { tmp=0.; goto fin; }
if( fabs(a) <= fabs(b)) tmp=a;
 else                   tmp=b;
fin:
return tmp;
}
/*--------------------------------------------------------------------*/
double riem_c(double c_l,double c_r,double spd)
{ 
double tmp;
if(spd >= 0.) tmp=c_l;
else          tmp=c_r;  
return tmp;
}
/*--------------------------------------------------------------------*/
