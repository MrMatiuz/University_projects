
/*CONTENTS
void inp_mesh(void);
void inp_gdf(void);
void derivs(void);
double step(void);
double flow_lr(void);
double flow_du(void);
void new_gdf(void);
*/



#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "func.h"

#define IN(i,j)   ( (i)*(JT+1)+(j) )
#define IC(i,j)   ( (i)*JT+(j) )
#define Ilr(i,j)  ( (i)*JT+(j) )
#define Idu(i,j)  ( (i)*(JT+1)+(j) )           


#define  GAM  1.4


#define max(a,b)  (((a) > (b)) ? (a) : (b))
#define min(a,b)  (((a) < (b)) ? (a) : (b))

double *gdf[6],*gdf_[6],*x,*y,*xc,*yc,*fllr[5],*fldu[5],
       *nlr[3],*tlr[2],*ndu[3],*tdu[2],*vol;
double *dxf[5],*dyf[5];
double *w[2];
double t,dt,t_max,coord,turb,Mach,Re_ref,u_ref,mu_ref;

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

// #define tau  mu_ref
#define tau  0.

double lxmin,lymin,mu_max,H_Z;
int n_stp,max_stp,inp_par,IT,JT,TW;

double RES;


void main(void)
{
 FILE *fpr;
 double d_t;

 fpr=fopen("par_i.dat","r");
 fscanf(fpr,"%d %d",&IT,&JT);
 fscanf(fpr,"%lg",&t_max);
 fscanf(fpr,"%d %d",&max_stp,&inp_par);
 fscanf(fpr,"%lg %lg",&Mach,&Re_ref);
 fscanf(fpr,"%lg %lg",&coord,&turb);
 fscanf(fpr,"%lg",&H_Z);
 fscanf(fpr,"%d",&TW);
 fclose(fpr);                                 

 t=0.;
 dt=1.e-7;

 inp_mesh();
 inp_gdf();
 n_stp=1;

 while( t < t_max )
   {
   double lcell,ki,bi;

   printf(" N_step=%d t=%e dt=%e RES=%e \n",n_stp,t,dt,RES);
   d_t=step();

   lcell=min(lxmin,lymin);
   ki=d_t;
   bi=0.25*lcell*lcell/mu_max;
   dt=min(bi,ki);
   dt=H_Z*dt;
   t+=dt;
   n_stp++;
   if(n_stp > max_stp) break;
   }
 out();
 out_tec();
}

//*****************************************************
double step(void)
{
 int i,j,k;
 double d_t,dt_lr,dt_du;

// BOUNDARY CONDITION:  ghost cells
// down:  no slip
   for(i=2; i < IT-2; i++)
    {
     gdf[0][IC(i,0)]=gdf[0][IC(i,1)]= TW/gdf[3][IC(i,2)];
     gdf[1][IC(i,0)]=gdf[1][IC(i,1)]=-gdf[1][IC(i,2)];
     gdf[2][IC(i,0)]=gdf[2][IC(i,1)]=-gdf[2][IC(i,2)];
     gdf[3][IC(i,0)]=gdf[3][IC(i,1)]= gdf[3][IC(i,2)];
     gdf[4][IC(i,0)]=gdf[4][IC(i,1)]= TW;      
    }

//up: in

// left  in

// right out
   for(j=0; j < JT; j++)
    {
    for(k=0; k < 5; k++)
     gdf[k][IC(IT-2,j)]=gdf[k][IC(IT-1,j)]=gdf[k][IC(IT-3,j)];
    }


 derivs();
 dt_lr=flow_lr();
 dt_du=flow_du();
 d_t=min(dt_lr,dt_du);
 new_gdf();
return d_t;
}

//*****************************************************

void inp_mesh(void)
{
 int i,j;
 double dx,dy;
 FILE *fpr;

 fpr=fopen("mesh.dat","r");
// inp coord
 fscanf(fpr,"%d %d",&IT,&JT);
 memo_on();

 for(i=2; i <= IT-2; i++)
  for(j=2; j <= JT-2; j++)
    fscanf(fpr,"%lg %lg",&x[IN(i,j)],&y[IN(i,j)]);

 lxmin=9.9e+9; lymin=9.9e+9;

// norms left_right
 for(i=2; i <= IT-2; i++)
  for(j=2; j < JT-2; j++)
   {
    dx=x[IN(i,j+1)]-x[IN(i,j)];
    dy=y[IN(i,j+1)]-y[IN(i,j)];
    nlr[0][Ilr(i,j)]=dy;
    nlr[1][Ilr(i,j)]=-dx;
    nlr[2][Ilr(i,j)]=sqrt(dx*dx+dy*dy);
    lxmin=min(lxmin,nlr[2][Ilr(i,j)]);
    tlr[0][Ilr(i,j)]=dx;
    tlr[1][Ilr(i,j)]=dy;
   }
// norms down_up
 for(i=2; i < IT-2; i++)
  for(j=2; j <= JT-2; j++)
   {
    dx=x[IN(i+1,j)]-x[IN(i,j)];
    dy=y[IN(i+1,j)]-y[IN(i,j)];
    ndu[0][Idu(i,j)]=-dy;
    ndu[1][Idu(i,j)]=dx;
    ndu[2][Idu(i,j)]=sqrt(dx*dx+dy*dy);
    lymin=min(lymin,ndu[2][Idu(i,j)]);
    tdu[0][Idu(i,j)]=dx;
    tdu[1][Idu(i,j)]=dy;
   }
// centers and volumes
 for(i=2; i < IT-2; i++)
  for(j=2; j < JT-2; j++)
   {
    dx=x[IN(i+1,j+1)]-x[IN(i,j)];
    dy=y[IN(i,j+1)]-y[IN(i+1,j)];
    vol[IC(i,j)]=0.5*dx*dy;
    dx=x[IN(i,j+1)]-x[IN(i+1,j)];
    dy=y[IN(i+1,j+1)]-y[IN(i,j)];
    vol[IC(i,j)]=vol[IC(i,j)]-0.5*dx*dy;
    xc[IC(i,j)]=0.25*(x[IN(i,j)]+x[IN(i+1,j)]+x[IN(i+1,j+1)]+x[IN(i,j+1)]);
    yc[IC(i,j)]=0.25*(y[IN(i,j)]+y[IN(i+1,j)]+y[IN(i+1,j+1)]+y[IN(i,j+1)]);
   }


// added 02.12.18
for(j=2; j < JT-2; j++)
  {
  xc[IC(IT-2,j)]=2.*xc[IC(IT-3,j)]-xc[IC(IT-4,j)];
  yc[IC(IT-2,j)]=2.*yc[IC(IT-3,j)]-yc[IC(IT-4,j)];
  xc[IC(IT-1,j)]=2.*xc[IC(IT-2,j)]-xc[IC(IT-3,j)];
  yc[IC(IT-1,j)]=2.*yc[IC(IT-2,j)]-yc[IC(IT-3,j)];

  xc[IC(1,j)]=2.*xc[IC(2,j)]-xc[IC(3,j)];
  yc[IC(1,j)]=2.*yc[IC(2,j)]-yc[IC(3,j)];
  xc[IC(0,j)]=2.*xc[IC(1,j)]-xc[IC(2,j)];
  yc[IC(0,j)]=2.*yc[IC(1,j)]-yc[IC(2,j)];
  }

for(i=2; i < IT-2; i++)
  {
  xc[IC(i,JT-2)]=2.*xc[IC(i,JT-3)]-xc[IC(i,JT-4)];
  yc[IC(i,JT-2)]=2.*yc[IC(i,JT-3)]-yc[IC(i,JT-4)];
  xc[IC(i,JT-1)]=2.*xc[IC(i,JT-2)]-xc[IC(i,JT-3)];
  yc[IC(i,JT-1)]=2.*yc[IC(i,JT-2)]-yc[IC(i,JT-3)];
                              
  xc[IC(i,1)]=2.*xc[IC(i,2)]-xc[IC(i,3)];
  yc[IC(i,1)]=2.*yc[IC(i,2)]-yc[IC(i,3)];
  xc[IC(i,0)]=2.*xc[IC(i,1)]-xc[IC(i,2)];
  yc[IC(i,0)]=2.*yc[IC(i,1)]-yc[IC(i,2)];
  }

 fclose(fpr);
}

//*****************************************************
void inp_gdf(void)
{
 int i,j,k;
 double d,u,v,p,e;
 double L_ref;
 FILE *fpr;
 d=1.; p=1.;  v=0.;
 u_ref=u=Mach*sqrt(GAM*p/d);
// !!!!!!!!!!!!!!!!
 L_ref=1.;
 mu_ref=u*L_ref*d/Re_ref;
 for(i=0; i < IT; i++)
  for(j=0; j < JT; j++)
   {
    k=IC(i,j);
    gdf[0][k]=d;
    gdf[1][k]=u;
    gdf[2][k]=v;
    gdf[3][k]=p;
    gdf[4][k]=p/d;
    gdf[5][k]=0.;
   }
 if(inp_par == 1)
  {
   fpr=fopen("reg_i.dat","rb");
   fread(&IT,sizeof(int),1,fpr);
   fread(&JT,sizeof(int),1,fpr);
   fread(&t,sizeof(double),1,fpr);
   fread(&dt,sizeof(double),1,fpr);
   for(i=2; i < IT-2; i++)
    for(j=2; j < JT-2; j++)
     for(k=0; k < 5; k++)
      {
      fread(&gdf[k][IC(i,j)],sizeof(double),1,fpr);
// !!!!!!!!!!!!!!!!Tempo 14118
      gdf_[k][IC(i,j)]=gdf[k][IC(i,j)]; 
      }
   fclose(fpr);
  }
// BOUNDARY CONDITION:  ghost cells

// down:  no slip
   for(i=2; i < IT-2; i++)
    {
     gdf[0][IC(i,0)]=gdf[0][IC(i,1)]= gdf[0][IC(i,2)];
     gdf[1][IC(i,0)]=gdf[1][IC(i,1)]=-gdf[1][IC(i,2)];
     gdf[2][IC(i,0)]=gdf[2][IC(i,1)]=-gdf[2][IC(i,2)];
     gdf[3][IC(i,0)]=gdf[3][IC(i,1)]= gdf[3][IC(i,2)];
     gdf[4][IC(i,0)]=gdf[4][IC(i,1)]= gdf[4][IC(i,2)];           //почему нет gdf[5]?
    }

//up: in

// left  in

// right out
   for(j=0; j < JT; j++)
    {
    for(k=0; k < 5; k++)
     gdf[k][IC(IT-2,j)]=gdf[k][IC(IT-1,j)]=gdf[k][IC(IT-3,j)];
    }

}

//*****************************************************
void derivs(void)
{
 int i,j,k,n;
 double f1,f2,f3,f4;


   for(i=2; i < IT-2; i++)
    for(j=2; j < JT-2; j++)
     {
     n=IC(i,j);
     for(k=0; k < 5; k++)
      {
      dxf[k][IC(i,j)]=0.;
      dyf[k][IC(i,j)]=0.;
      f1=0.25*(gdf[k][IC(i-1,j-1)]+gdf[k][IC(i,j-1)]+
               gdf[k][IC(i,j)]+gdf[k][IC(i-1,j)]);             //???
      f2=0.25*(gdf[k][IC(i,j-1)]+gdf[k][IC(i+1,j-1)]+
               gdf[k][IC(i+1,j)]+gdf[k][IC(i,j)]);
      f3=0.25*(gdf[k][IC(i,j)]+gdf[k][IC(i+1,j)]+
               gdf[k][IC(i+1,j+1)]+gdf[k][IC(i,j+1)]);
      f4=0.25*(gdf[k][IC(i-1,j)]+gdf[k][IC(i,j)]+
               gdf[k][IC(i,j+1)]+gdf[k][IC(i-1,j+1)]);

      dxf[k][IC(i,j)]-=0.5*(f1+f2)*ndu[0][Idu(i,j)];           //???
      dyf[k][IC(i,j)]-=0.5*(f1+f2)*ndu[1][Idu(i,j)];

      dxf[k][IC(i,j)]+=0.5*(f2+f3)*nlr[0][Ilr(i+1,j)];
      dyf[k][IC(i,j)]+=0.5*(f2+f3)*nlr[1][Ilr(i+1,j)];

      dxf[k][IC(i,j)]+=0.5*(f3+f4)*ndu[0][Idu(i,j+1)];
      dyf[k][IC(i,j)]+=0.5*(f3+f4)*ndu[1][Idu(i,j+1)];


      dxf[k][IC(i,j)]-=0.5*(f1+f4)*nlr[0][Ilr(i,j)];
      dyf[k][IC(i,j)]-=0.5*(f1+f4)*nlr[1][Ilr(i,j)];

      dxf[k][IC(i,j)]=dxf[k][IC(i,j)]/vol[IC(i,j)];
      dyf[k][IC(i,j)]=dyf[k][IC(i,j)]/vol[IC(i,j)];
      }

      w[0][IC(i,j)]=0.5*tau*( gdf[1][IC(i,j)]*gdf[1][IC(i,j)]*dxf[0][IC(i,j)]+
                            2*gdf[0][IC(i,j)]*gdf[1][IC(i,j)]*dxf[1][IC(i,j)]+
                            dxf[3][IC(i,j)] + gdf[1][IC(i,j)]*gdf[2][IC(i,j)]*dyf[0][IC(i,j)]+
                            gdf[0][IC(i,j)]*gdf[2][IC(i,j)]*dyf[1][IC(i,j)]+
                            gdf[0][IC(i,j)]*gdf[1][IC(i,j)]*dyf[2][IC(i,j)])/gdf[0][IC(i,j)];

      w[1][IC(i,j)]=0.5*tau*( gdf[1][IC(i,j)]*gdf[2][IC(i,j)]*dxf[0][IC(i,j)]+
                            gdf[0][IC(i,j)]*gdf[2][IC(i,j)]*dxf[1][IC(i,j)]+
                            gdf[0][IC(i,j)]*gdf[1][IC(i,j)]*dxf[2][IC(i,j)]+
                            gdf[2][IC(i,j)]*gdf[2][IC(i,j)]*dyf[0][IC(i,j)]+
                            2*gdf[0][IC(i,j)]*gdf[2][IC(i,j)]*dyf[2][IC(i,j)]+
                            dyf[3][IC(i,j)] )/gdf[0][IC(i,j)];
     }

}
//*****************************************************

double flow_lr(void)
{
 int i,j,k,l,lom,alarm;

 double dx,dy,dksi1,dksi2,delf1,delf2,delf;
 double db,ub,vb,pb,e,NV1,TV1,NV2,TV2;
 double w1, w2;
 double t11,t12,t21,t22,q1,q2;
 double mum,mut,muef,kapm,kapt,kapef;
 double ds,ff[5],dxflr[5],dyflr[5];
 double dt1,dt2,dt_max,dt_lr;
 static double par[32],rqp[3][2],r_o,u_o,v_o,p_o,om,s1,s2,s3;

 lom=0;  om=0.; alarm=0;
 mu_max=0.;
 dt_lr=1.e17;

 for(i=2; i <= IT-2; i++)
  for(j=2; j < JT-2; j++)
   {
// left  values for riem
	dx = xc[IC(i-1,j)] - xc[IC(i - 2,j)];
	dy = yc[IC(i-1,j)] - yc[IC(i - 2,j)];
	dksi1 = sqrt(dx*dx+dy*dy);
	dx = xc[IC(i,j)] - xc[IC(i - 1,j)];
	dy = yc[IC(i,j)] - yc[IC(i - 1,j)];
	dksi2 = sqrt(dx*dx+dy*dy);
	for( k=0; k < 5; k++)
	{
		delf1 = (gdf[k][IC(i-1,j)] - gdf[k][IC(i-2,j)]) / dksi1;
		delf2 = (gdf[k][IC(i,j)] - gdf[k][IC(i-1,j)]) / dksi2;
		delf = min_mod(delf1,delf2);

//if(i==2 || i== IT-2) delf = 0.;

		ff[k]=gdf[k][IC(i-1,j)]+0.5*dksi2*delf;
	}
     NV1=(nlr[0][Ilr(i,j)]*ff[1]+nlr[1][Ilr(i,j)]*ff[2])/nlr[2][Ilr(i,j)];
     TV1=(tlr[0][Ilr(i,j)]*ff[1]+tlr[1][Ilr(i,j)]*ff[2])/nlr[2][Ilr(i,j)];
     rqp[0][0]=ff[0];
     rqp[1][0]=NV1;
     rqp[2][0]=ff[3];

// right values for riem
        dksi1=dksi2;
	dx = xc[IC(i+1,j)] - xc[IC(i,j)];
	dy = yc[IC(i+1,j)] - yc[IC(i,j)];
	dksi2 = sqrt(dx*dx+dy*dy);
	for( k=0; k < 5; k++)
	{
		delf1=(gdf[k][IC(i,j)]-gdf[k][IC(i-1,j)])/dksi1;
		delf2=(gdf[k][IC(i+1,j)]-gdf[k][IC(i,j)])/dksi2;
		delf=min_mod(delf1,delf2);

//if(i==2 || i== IT-2) delf = 0.;

		ff[k]=gdf[k][IC(i,j)]-0.5*dksi1*delf;
	}
     NV2=(nlr[0][Ilr(i,j)]*ff[1]+nlr[1][Ilr(i,j)]*ff[2])/nlr[2][Ilr(i,j)];
     TV2=(tlr[0][Ilr(i,j)]*ff[1]+tlr[1][Ilr(i,j)]*ff[2])/nlr[2][Ilr(i,j)];
     rqp[0][1]=ff[0];
     rqp[1][1]=NV2;
     rqp[2][1]=ff[3];

       for(l=0; l < 2; l++)
        for(k=0; k < 3; k++)
         par[k+3*l]=rqp[k][l];

       riemi(par,&alarm);
        if(alarm )
          {
          printf("ALARM_RIEM FROM l_r : i=%d j=%d  \n",i,j);
          printf("L: %e %e %e \n",rqp[0][0],rqp[1][0],rqp[2][0]);
          printf("R: %e %e %e \n",rqp[0][1],rqp[1][1],rqp[2][1]);
          exit(-1);
          }

     r_o=par[6];  u_o=par[7] ;   p_o=par[8] ;
     s1 =par[9];  s2 =par[10];   s3 =par[11];

     v_o=TV1;
     if(s2 < om) v_o=TV2;

     db=r_o; pb=p_o;
     ub=(u_o*nlr[0][Ilr(i,j)]+v_o*tlr[0][Ilr(i,j)])/nlr[2][Ilr(i,j)];
     vb=(u_o*nlr[1][Ilr(i,j)]+v_o*tlr[1][Ilr(i,j)])/nlr[2][Ilr(i,j)];
     e=pb/(GAM-1.)+0.5*db*(ub*ub+vb*vb);

// Eu_flux
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!231118
     if(s2 < om)
     {     
     w1 = w[0][IC(i,j)];
     w2 = w[1][IC(i,j)];
     } 
     else
     {     
     w1 = w[0][IC(i-1,j)];
     w2 = w[1][IC(i-1,j)];
     } 
if(i == 2)
	{
	w1 = w[0][IC(i,j)];
	w2 = w[1][IC(i,j)];
	}
if(i == IT-2)
	{
	w1 = w[0][IC(i-1,j)];
	w2 = w[1][IC(i-1,j)];
	}

     fllr[0][Ilr(i,j)]=db*(ub - w1)*nlr[0][Ilr(i,j)]+db*(vb-w2)*nlr[1][Ilr(i,j)];
     fllr[1][Ilr(i,j)]=(db*ub*(ub - w1)+pb)*nlr[0][Ilr(i,j)]+db*(ub - w1)*vb*nlr[1][Ilr(i,j)];
     fllr[2][Ilr(i,j)]=db*ub*(vb-w2)*nlr[0][Ilr(i,j)]+(db*vb*(vb-w2)+pb)*nlr[1][Ilr(i,j)];
     fllr[3][Ilr(i,j)]=(ub - w1)*(e+pb)*nlr[0][Ilr(i,j)]+(vb-w2)*(e+pb)*nlr[1][Ilr(i,j)];

// NS_flux
// !!!!!!!!!!!!!!!!!!!!
     muef=mu_ref*pow((pb/db),0.76);
     kapef=muef/0.72;
     kapef*=3.5;
     if(mu_max < muef) mu_max=muef;

     for(k=0; k < 5; k++)
       {
        dxflr[k]=0.5*(dxf[k][IC(i-1,j)]+dxf[k][IC(i,j)]);
        dyflr[k]=0.5*(dyf[k][IC(i-1,j)]+dyf[k][IC(i,j)]);
if(i == 2)
	{
	dxflr[k] = dxf[k][IC(i, j)];
	dyflr[k] = dyf[k][IC(i, j)];
	}
if(i == IT-2)
	{
	dxflr[k] = dxf[k][IC(i - 1, j)];
	dyflr[k] = dyf[k][IC(i - 1, j)];
	}
       }

     t11=2.*muef*(2.*dxflr[1]-dyflr[2])/3.;
     t12=muef*(dyflr[1]+dxflr[2]);
     t21=t12;
     t22=2.*muef*(2.*dyflr[2]-dxflr[1])/3.;
     q1=kapef*dxflr[4];
     q2=kapef*dyflr[4];


     fllr[1][Ilr(i,j)]-=t11*nlr[0][Ilr(i,j)]+t12*nlr[1][Ilr(i,j)];
     fllr[2][Ilr(i,j)]-=t21*nlr[0][Ilr(i,j)]+t22*nlr[1][Ilr(i,j)];
     fllr[3][Ilr(i,j)]-=(ub*t11+vb*t12+q1)*nlr[0][Ilr(i,j)]+
                        (ub*t21+vb*t22+q2)*nlr[1][Ilr(i,j)];

// calc time-step
if(i < IT-2)
{
     dt1=1.e17;
     if( s1 < 0.)
     dt1=dksi1/fabs(s1);
     dt2=1.e17;
     if( s3 > 0.)
     dt2=dksi2/fabs(s3);
     dt_max=min(dt1,dt2);
     dt_lr=min(dt_max,dt_lr);
}

   }
return dt_lr;
}

//*****************************************************
double min_mod(double a, double b)
{
 double tmp;
 tmp=0.;

 if( a*b > 0.)
  {
   tmp=a;
   if( fabs(a) > fabs(b) )  tmp=b;
  }

 return tmp;
}

//*****************************************************

double flow_du(void)
{
 int i,j,k,l,lom,alarm;
 double dx,dy,deta1,deta2,delf1,delf2,delf;
 double db,ub,vb,pb,e,NV1,TV1,NV2,TV2;
 double w1, w2;
 double t11,t12,t21,t22,q1,q2;
 double mum,mut,muef,kapm,kapt,kapef;
 double ds,ff[5],dxfdu[5],dyfdu[5];
 double dt1,dt2,dt_max,dt_du;
 static double par[32],rqp[3][2],r_o,u_o,v_o,p_o,om,s1,s2,s3;

 lom=0;  om=0.;
 mu_max=0.;
 dt_du=1.e17;
 for(i=2; i < IT-2; i++)
  for(j=2; j <= JT-2; j++)
   {
// left(down)  values for riem
	dx=xc[IC(i,j-1)]-xc[IC(i,j-2)];
	dy=yc[IC(i,j-1)]-yc[IC(i,j-2)];
	deta1=sqrt(dx*dx+dy*dy);
	dx=xc[IC(i,j)]-xc[IC(i,j-1)];
	dy=yc[IC(i,j)]-yc[IC(i,j-1)];
	deta2=sqrt(dx*dx+dy*dy);
	for( k=0; k < 5; k++)
	{
		delf1=(gdf[k][IC(i,j-1)]-gdf[k][IC(i,j-2)])/deta1;
		delf2=(gdf[k][IC(i,j)]-gdf[k][IC(i,j-1)])/deta2;
		delf=min_mod(delf1,delf2);

//if(j==2 || j== JT-2) delf=0.;

		ff[k]=gdf[k][IC(i,j-1)]+0.5*deta2*delf;
	}
     NV1=(ndu[0][Idu(i,j)]*ff[1]+ndu[1][Idu(i,j)]*ff[2])/ndu[2][Idu(i,j)];
     TV1=(tdu[0][Idu(i,j)]*ff[1]+tdu[1][Idu(i,j)]*ff[2])/ndu[2][Idu(i,j)];
     rqp[0][0]=ff[0];
     rqp[1][0]=NV1;
     rqp[2][0]=ff[3];
// right (up) values for riem
	deta1=deta2;
	dx=xc[IC(i,j+1)]-xc[IC(i,j)];
	dy=yc[IC(i,j+1)]-yc[IC(i,j)];
	deta2=sqrt(dx*dx+dy*dy);
	for( k=0; k < 4; k++)
	{
		delf1=(gdf[k][IC(i,j)]-gdf[k][IC(i,j-1)])/deta1;
		delf2=(gdf[k][IC(i,j+1)]-gdf[k][IC(i,j)])/deta2;
		delf=min_mod(delf1,delf2);

//if(j==2 || j==JT-2) delf=0.;

		ff[k]=gdf[k][IC(i,j)]-0.5*deta1*delf;
	}
     NV2=(ndu[0][Idu(i,j)]*ff[1]+ndu[1][Idu(i,j)]*ff[2])/ndu[2][Idu(i,j)];
     TV2=(tdu[0][Idu(i,j)]*ff[1]+tdu[1][Idu(i,j)]*ff[2])/ndu[2][Idu(i,j)];
     rqp[0][1]=ff[0];
     rqp[1][1]=NV2;
     rqp[2][1]=ff[3];

// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     if(j==2)
      {
      NV1=-NV2;
      TV1=-TV2;
      rqp[0][0]=rqp[0][1];
      rqp[1][0]=-rqp[1][1];
      rqp[2][0]=rqp[2][1];
      }


       for(l=0; l < 2; l++)
        for(k=0; k < 3; k++)
         par[k+3*l]=rqp[k][l];

       riemi(par,&alarm);
        if(alarm )
          {
          printf("ALARM_RIEM FROM d_u : i=%d j=%d  \n",i,j);
          printf("L: %e %e %e \n",rqp[0][0],rqp[1][0],rqp[2][0]);
          printf("R: %e %e %e \n",rqp[0][1],rqp[1][1],rqp[2][1]);
          exit(-1);
          }

     r_o=par[6];  u_o=par[7] ;   p_o=par[8] ;
     s1 =par[9];  s2 =par[10];   s3 =par[11];

     v_o=TV1;
     if(s2  < om) v_o=TV2;
     db=r_o; pb=p_o;
     ub=(u_o*ndu[0][Idu(i,j)]+v_o*tdu[0][Idu(i,j)])/ndu[2][Idu(i,j)];
     vb=(u_o*ndu[1][Idu(i,j)]+v_o*tdu[1][Idu(i,j)])/ndu[2][Idu(i,j)];
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     if(j==2) ub=vb=0.;

     e=pb/(GAM-1.)+0.5*db*(ub*ub+vb*vb);


// Eu_flux
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!231118
     if(s2 < om)
     {     
     w1 = w[0][IC(i,j)];
     w2 = w[1][IC(i,j)];
     } 
     else
     {     
     w1 = w[0][IC(i,j-1)];
     w2 = w[1][IC(i,j-1)];
     } 
//     w1 = 0.5*(w[0][IC(i,j-1)]+w[0][IC(i,j)]);
//     w2 = 0.5*(w[1][IC(i,j-1)]+w[1][IC(i,j)]);
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!141118
if(j==2)
 {
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!231118
// w1=w[0][IC(i,j)];
// w2=w[1][IC(i,j)];
w1=0.;
w2=0.;
 }
if(j==JT-2)
 {
 w1=w[0][IC(i,j-1)];
 w2=w[1][IC(i,j-1)];
 }


     fldu[0][Idu(i,j)]=db*(ub-w1)*ndu[0][Idu(i,j)]+db*(vb-w2)*ndu[1][Idu(i,j)];
     fldu[1][Idu(i,j)]=(db*ub*(ub-w1)+pb)*ndu[0][Idu(i,j)]+db*(ub-w1)*vb*ndu[1][Idu(i,j)];
     fldu[2][Idu(i,j)]=db*ub*(vb-w2)*ndu[0][Idu(i,j)]+(db*vb*(vb-w2)+pb)*ndu[1][Idu(i,j)];
     fldu[3][Idu(i,j)]=(ub-w1)*(e+pb)*ndu[0][Idu(i,j)]+(vb-w2)*(e+pb)*ndu[1][Idu(i,j)];

// NS_flux
     muef=mu_ref*pow((pb/db),0.76);
     kapef=muef/0.72;
     kapef*=3.5;
     if(mu_max < muef) mu_max=muef;
     for(k=0; k < 5; k++)
       {
        dxfdu[k]=0.5*(dxf[k][IC(i,j-1)]+dxf[k][IC(i,j)]);
        dyfdu[k]=0.5*(dyf[k][IC(i,j-1)]+dyf[k][IC(i,j)]);
if(j==2)
 {
 dxfdu[k]=dxf[k][IC(i,j)];
 dyfdu[k]=dyf[k][IC(i,j)];
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!for B/L ONLY
// ??????????????????????????????????
if(i > 2)
 {
  dxfdu[1]=0.;
  dxfdu[2]=0.;
 }
 }
if(j==JT-2)
 {
 dxfdu[k]=dxf[k][IC(i,j-1)];
 dyfdu[k]=dyf[k][IC(i,j-1)];
 }
       }

     t11=2.*muef*(2.*dxfdu[1]-dyfdu[2])/3.;
     t12=muef*(dyfdu[1]+dxfdu[2]);
     t21=t12;
     t22=2.*muef*(2.*dyfdu[2]-dxfdu[1])/3.;
     q1=kapef*dxfdu[4];
     q2=kapef*dyfdu[4];

     fldu[1][Idu(i,j)]-=t11*ndu[0][Idu(i,j)]+t12*ndu[1][Idu(i,j)];
     fldu[2][Idu(i,j)]-=t21*ndu[0][Idu(i,j)]+t22*ndu[1][Idu(i,j)];
     fldu[3][Idu(i,j)]-=(ub*t11+vb*t12+q1)*ndu[0][Idu(i,j)]+
                        (ub*t21+vb*t22+q2)*ndu[1][Idu(i,j)];

// calc time-step
if( j < JT-2)
{
     dt1=1.e17;
     if(s1 < 0.)
     dt1=deta1/fabs(s1);
     dt2=1.e17;
     if(s3 > 0.)
     dt2=deta2/fabs(s3);
     dt_max=min(dt1,dt2);
     dt_du=min(dt_max,dt_du);
}

   }
return dt_du;
}


//*****************************************************

void new_gdf(void)
{
 int i,j,n,k;
 double p,e;
 static double rhs[5],a0[5],a1[5],a_[5],fl[5];

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 rhs[0]=rhs[1]=rhs[2]=rhs[3]=rhs[4]=0.;
 RES=0.;

 for(i=2; i < IT-2; i++)
  for(j=2; j < JT-2; j++)
   {
    n=IC(i,j);
//!!!!!!!!!!!!!!!!!!!!!!!!!! 221118
//    tau=mu_ref*pow((gdf[3][n]/gdf[0][n]),0.76)/gdf[3][n];
    a0[0]=gdf[0][n];
    a0[1]=gdf[0][n]*gdf[1][n];
    a0[2]=gdf[0][n]*gdf[2][n];
    e=gdf[3][n]/(GAM-1.)+0.5*gdf[0][n]*
      (gdf[1][n]*gdf[1][n]+gdf[2][n]*gdf[2][n]);
    a0[3]=e;

    a_[0]=gdf_[0][n];
    a_[1]=gdf_[0][n]*gdf_[1][n];
    a_[2]=gdf_[0][n]*gdf_[2][n];
    e=gdf_[3][n]/(GAM-1.)+0.5*gdf_[0][n]*
      (gdf_[1][n]*gdf_[1][n]+gdf_[2][n]*gdf_[2][n]);
    a_[3]=e;

    for(k=0; k < 4; k++)
     {
     fl[k]=(fllr[k][Ilr(i+1,j)]-fllr[k][Ilr(i,j)]+fldu[k][Idu(i,j+1)]-fldu[k][Idu(i,j)]);
     a1[k]=a0[k]-dt*fl[k]/vol[n]+dt*rhs[k]-tau*(a_[k]-2.*a0[k])/dt;
     a1[k]=a1[k]/(1.+tau/dt);
     gdf_[k][n]=gdf[k][n];
     }
    RES=max(RES,fabs(fl[0]));
    gdf[0][n]=a1[0];
    gdf[1][n]=a1[1]/a1[0];
    gdf[2][n]=a1[2]/a1[0];
    p=(GAM-1.)*(a1[3]-0.5*(a1[1]*gdf[1][n]+a1[2]*gdf[2][n]));
    gdf[3][n]=p;
    gdf[4][n]=p/gdf[0][n];
   }
}

//*****************************************************
/*
double min(double m1,double m2)
{
    double res;
    if (m1<=m2) 
        res = m1;
    else
        res = m2;
    return res;
}
double max(double m1,double m2)
{
    if (m1>=m2) 
        res = m1;
    else
        res = m2;
    return res;
}
*/

