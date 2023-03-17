#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "func.h"

#define IN(i,j)   ( (i)*(JT+1)+(j) )
#define IC(i,j)   ( (i)*JT+(j) )
#define Ilr(i,j)  ( (i)*JT+(j) )
#define Idu(i,j)  ( (i)*(JT+1)+(j) )           

#define  GAM  1.4


extern double *gdf[6],*gdf_[6],*x,*y,*xc,*yc,*fllr[5],*fldu[5],
       *nlr[3],*tlr[2],*ndu[3],*tdu[2],*vol;
extern double *dxf[5],*dyf[5];
extern double *w[2];
extern double t,dt,t_max,coord,turb,Mach,Re_ref,u_ref,mu_ref;

extern int n_out, i_app,i_flow,n_stp,max_stp,inp_par,IT,JT;
extern double tau;

#define max(a,b)  (((a) > (b)) ? (a) : (b))
#define min(a,b)  (((a) < (b)) ? (a) : (b))

#define bet 2.


void fllrRusanov(int i, int j, double *frul, double *frur,double &db, double &ub, double &vb, double &pb,double *data)
{
	double ds,e,NV1,NV2, ff[5];
        double w1,w2;
	int i_lr = Ilr(i,j);
        ds=nlr[2][i_lr];
	for(int k=0; k < 4; k++)
		ff[k]=frul[k]; 

	NV1=(nlr[0][i_lr]*ff[1]+nlr[1][i_lr]*ff[2])/ds;
	double fl[5],fr[5],cl,cr,lam;
	db=ff[0];  ub=ff[1]; vb=ff[2];  pb=ff[3];
	e=pb/(GAM-1.)+0.5*db*(ub*ub+vb*vb);



     w1 = 0.5*(w[0][IC(i-1,j)]+w[0][IC(i,j)]);
     w2 = 0.5*(w[1][IC(i-1,j)]+w[1][IC(i,j)]);

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

    fl[0] =( db*(ub - w1)*nlr[0][Ilr(i,j)]+db*(vb-w2)*nlr[1][Ilr(i,j)] )/ds;
    fl[1] =( (db*ub*(ub - w1)+pb)*nlr[0][Ilr(i,j)]+db*(ub - w1)*vb*nlr[1][Ilr(i,j)] )/ds;
    fl[2] =( db*ub*(vb-w2)*nlr[0][Ilr(i,j)]+(db*vb*(vb-w2)+pb)*nlr[1][Ilr(i,j)] )/ds;
    fl[3] =( (ub - w1)*(e+pb)*nlr[0][Ilr(i,j)]+(vb-w2)*(e+pb)*nlr[1][Ilr(i,j)] )/ds;

	cl=sqrt(GAM*pb/db); 
	double ql[] = {db, db * ub, db * vb, e};


	for(int k=0; k < 4; k++)
		ff[k]=frur[k]; 

	db=ff[0];  ub=ff[1]; vb=ff[2];  pb=ff[3];
	NV2=(nlr[0][i_lr]*ff[1]+nlr[1][i_lr]*ff[2])/ds;
	e=pb/(GAM-1.)+0.5*db*(ub*ub+vb*vb);

	fr[0]=( db*(ub - w1)*nlr[0][Ilr(i,j)]+db*(vb-w2)*nlr[1][Ilr(i,j)] )/ds;              
	fr[1]=( (db*ub*(ub - w1)+pb)*nlr[0][Ilr(i,j)]+db*(ub - w1)*vb*nlr[1][Ilr(i,j)] )/ds; 
	fr[2]=( db*ub*(vb-w2)*nlr[0][Ilr(i,j)]+(db*vb*(vb-w2)+pb)*nlr[1][Ilr(i,j)] )/ds;     
	fr[3]=( (ub - w1)*(e+pb)*nlr[0][Ilr(i,j)]+(vb-w2)*(e+pb)*nlr[1][Ilr(i,j)] )/ds;      

	cr=sqrt(GAM*pb/db); 
	double qr[] = {db, db * ub, db * vb, e};
	lam=max( fabs(NV1-cl),fabs(NV2+cr) );
        data[9]=NV1-cl;
        data[11]=NV2+cr;
	for(int k=0; k < 4; k++) 
	  fllr[k][i_lr]=(0.5*(fl[k]+fr[k])-0.5*bet*lam*(qr[k]-ql[k]))*ds;

	db = 0.5 * (gdf[0][IC(i-1,j)] + gdf[0][IC(i,j)]);
	ub = 0.5 * (gdf[1][IC(i-1,j)] + gdf[1][IC(i,j)]);
	vb = 0.5 * (gdf[2][IC(i-1,j)] + gdf[2][IC(i,j)]);
	pb = 0.5 * (gdf[3][IC(i-1,j)] + gdf[3][IC(i,j)]);
}


void flduRusanov(int i, int j, double *frul, double *frur, double &db, double &ub, double &vb, double &pb, double *data)
{
	double ds,e,NV1,NV2, ff[5];
	int i_du = Idu(i, j);
        ds=ndu[2][i_du];
       	for(int k=0; k < 4; k++)
	 ff[k]=frul[k]; 

	double fl[5],fr[5],cl,cr,lam;
	db=ff[0];  ub=ff[1]; vb=ff[2];  pb=ff[3];
	NV1=(ndu[0][i_du]*ff[1]+ndu[1][i_du]*ff[2])/ds;
	e=pb/(GAM-1.)+0.5*db*(ub*ub+vb*vb);
	fl[0]=(db*ub*ndu[0][i_du]+db*vb*ndu[1][i_du])/ds;
	fl[1]=((db*ub*ub+pb)*ndu[0][i_du]+db*ub*vb*ndu[1][i_du])/ds;
	fl[2]=(db*ub*vb*ndu[0][i_du]+(db*vb*vb+pb)*ndu[1][i_du])/ds;
	fl[3]=(ub*(e+pb)*ndu[0][i_du]+vb*(e+pb)*ndu[1][i_du])/ds;
	cl=sqrt(GAM*pb/db); 
	double ql[] = {db, db * ub, db * vb, e};

	for(int k=0; k < 4; k++)
	ff[k]=frur[k]; 

	db=ff[0];  ub=ff[1]; vb=ff[2];  pb=ff[3];
	NV2=(ndu[0][i_du]*ff[1]+ndu[1][i_du]*ff[2])/ds;
	e=pb/(GAM-1.)+0.5*db*(ub*ub+vb*vb);
	fr[0]=(db*ub*ndu[0][i_du]+db*vb*ndu[1][i_du])/ds;
	fr[1]=((db*ub*ub+pb)*ndu[0][i_du]+db*ub*vb*ndu[1][i_du])/ds;
	fr[2]=(db*ub*vb*ndu[0][i_du]+(db*vb*vb+pb)*ndu[1][i_du])/ds;
	fr[3]=(ub*(e+pb)*ndu[0][i_du]+vb*(e+pb)*ndu[1][i_du])/ds;
	cr=sqrt(GAM*pb/db); 
	double qr[] = {db, db * ub, db * vb, e};
	lam=max( fabs(NV1-cl),fabs(NV2+cr) );
        data[9]=NV1-cl;
        data[11]=NV2+cr;
	for(int k=0; k < 4; k++) 
         { 
	  fldu[k][i_du]=(0.5*(fl[k]+fr[k])-0.5*bet*lam*(qr[k]-ql[k]))*ds;
         }

	db = 0.5 * (gdf[0][IC(i,j-1)] + gdf[0][IC(i,j)]);
	ub = 0.5 * (gdf[1][IC(i,j-1)] + gdf[1][IC(i,j)]);
        vb = 0.5 * (gdf[2][IC(i,j-1)] + gdf[2][IC(i,j)]);
	pb = 0.5 * (gdf[3][IC(i,j-1)] + gdf[3][IC(i,j)]);
}
