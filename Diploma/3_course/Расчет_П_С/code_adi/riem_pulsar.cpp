
//  Солверы для расчета задачи о распаде разрыва
//  Текущий: неявный метод Годунова

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "riem.h"

void riemi(double* rqp_io, int *alarm)
{

	double om = 0;

	double c[2],rc[2],du,ff,df;
	double ss[2][2],rr[2],cc[2];
	//definition for ff_df
	double rc_1;
	double f[2];

	int it;
	double e_pmin,u_vc,p0,p1,pp,uu,ai,tmp1,tmp2;
        
        int lom=0;
        *alarm=0;

	// i==0 //
	{
		c[0]=sqrt(g*rqp_io[2]/rqp_io[0]);
		rc[0]=rqp_io[0]*c[0];
	}
	// i==1 //
	{
		c[1]=sqrt(g*rqp_io[5]/rqp_io[3]);
		rc[1]=rqp_io[3]*c[1];
	}
	e_pmin=min(rqp_io[2],rqp_io[5])*ee;

	/* check existance */
	u_vc=-c[0]*gc-c[1]*gc;
	du=rqp_io[1]-rqp_io[4];
	if(du <= u_vc)
	{
	printf("\n vacuem - no solution \n");
	printf(" Left  %e  %e  %e  \n",rqp_io[0],rqp_io[1],rqp_io[2]);
	printf(" Right %e  %e  %e  \n",rqp_io[3],rqp_io[4],rqp_io[5]);
        *alarm=77;
        return;
//	fflush(stdout);
//	exit (-1);
	}

	/* inicial guess */
	p1=(rqp_io[2]*rc[1]+rqp_io[5]*rc[0]+du*rc[0]*rc[1])/(rc[0]+rc[1]);
	p1=max(p1,e_pmin);
	it=1;

	/* Newton iter */
	do
	{
		if( it > maxit )
		{
		printf(" ITER > MAXIT \n");  
        	printf(" Left  %e  %e  %e  \n",rqp_io[0],rqp_io[1],rqp_io[2]);
	        printf(" Right %e  %e  %e  \n",rqp_io[3],rqp_io[4],rqp_io[5]);
                *alarm=78;
                return;
		}
		p0=p1;
		//ff_df(p0)
		{
			df=0.;
			// i==0 //
			{
				pp=p0/rqp_io[2];
				if(pp > 1.)
				{
					rc_1=1.f/rc[0];
					tmp1=1.f/sqrt(ga*pp+gb);
					f[0]=(p0-rqp_io[2])*tmp1*rc_1;
					tmp2=(g+1.)*pp+3.*g-1.;
					df=df+0.25*tmp2*(tmp1*tmp1*tmp1)*rc_1*gg;
				}
				else
				{
					tmp1=pow(pp,gb);
					f[0]=gc*c[0]*(tmp1-1.);
					df=df+c[0]*tmp1*gg/p0;
				}
			}
			// i==1 //
			{
				pp=p0/rqp_io[5];
				if(pp > 1.)
				{
					rc_1=1./rc[1];
					tmp1=1./sqrt(ga*pp+gb);
					f[1]=(p0-rqp_io[5])*tmp1*rc_1;
					tmp2=(g+1.)*pp+3.*g-1.;
					df=df+0.25*tmp2*(tmp1*tmp1*tmp1)*rc_1*gg;
				}
				else
				{
					tmp1=pow(pp,gb);
					f[1]=gc*c[1]*(tmp1-1.);
					df=df+c[1]*tmp1*gg/p0;
				}
			}
			ff=f[0]+f[1]-du;
			rqp_io[10]=rqp_io[1]-f[0];
		}
		p1=p0-ff/df;
		it++;
	}
	while ( fabs(p1-p0) > e_pmin );
 
	/* wave coeff. and dens */
	pp=p0; uu=rqp_io[10];
	// i==0 //
	{
		ai=2.*( (double) 0.0 -0.5);
		if(p0 < rqp_io[2]*estr)
		{
			ss[0][0]=rqp_io[1]+ai*c[0];
			cc[0]=c[0]-0.5*ai*(g-1.)*(rqp_io[1]-uu);
			ss[1][0]=uu+ai*cc[0];
			rr[0]=g*pp/(cc[0]*cc[0]);
		}
		else
		{
			tmp1=(g+1.)*pp+(g-1.)*rqp_io[2];
			tmp2=(g-1.)*pp+(g+1.)*rqp_io[2];
			rr[0]=rqp_io[0]*tmp1/tmp2;
			tmp1=rqp_io[0]-rr[0];
			tmp2=rqp_io[0]*rqp_io[1]-rr[0]*uu;
			ss[0][0]=tmp2/tmp1;
			ss[1][0]=ss[0][0];
		}
	}
	// i==1 //
	{
		ai=2.*( (double) 1.0 -0.5);
		if(p0 < rqp_io[5]*estr)
		{
			ss[0][1]=rqp_io[4]+ai*c[1];
			cc[1]=c[1]-0.5*ai*(g-1.)*(rqp_io[4]-uu);
			ss[1][1]=uu+ai*cc[1];
			rr[1]=g*pp/(cc[1]*cc[1]);
		}
		else
		{
			tmp1=(g+1.)*pp+(g-1.)*rqp_io[5];
			tmp2=(g-1.)*pp+(g+1.)*rqp_io[5];
			rr[1]=rqp_io[3]*tmp1/tmp2;
			tmp1=rqp_io[3]-rr[1];
			tmp2=rqp_io[3]*rqp_io[4]-rr[1]*uu;
			ss[0][1]=tmp2/tmp1;
			ss[1][1]=ss[0][1];
		}
	}
 
	/* data on line */
	if(lom == 1) om=ss[0][0];
	if(lom == 2) om=uu;
	if(lom == 3) om=ss[0][1];
	rqp_io[9]=ss[0][0];
	rqp_io[11]=ss[0][1];

	/* 1 */
	if( om <= ss[0][0] )
	{
		rqp_io[6]=rqp_io[0]; rqp_io[7]=rqp_io[1]; rqp_io[8]=rqp_io[2]; goto rtn_o; 
	}

	/* 2 */
	if( ( om > ss[0][0] ) && ( om <= ss[1][0]) )
	{
		tmp1=1./(g+1.);
		cc[0]=(2.*c[0]+(g-1.)*(rqp_io[1]-om))*tmp1;
		rqp_io[7]=om+cc[0];
		rqp_io[6]=rqp_io[0]*pow( (cc[0]/c[0]),gc );
		rqp_io[8]=cc[0]*cc[0]*rqp_io[6]/g;
	goto rtn_o;
	}

	/* 3 */
	if( (om > ss[1][0]) && (om <= rqp_io[10]) ){ rqp_io[6]=rr[0]; rqp_io[7]=uu; rqp_io[8]=pp; goto rtn_o; }

	/* 4 */
	if( (om > rqp_io[10]) && (om <= ss[1][1]) ){ rqp_io[6]=rr[1]; rqp_io[7]=uu; rqp_io[8]=pp; goto rtn_o; }

	/* 5 */
	if( (om > ss[1][1]) && (om <= ss[0][1]) )
	{
		tmp2=1./(g+1.);
		cc[1]=(2.*c[1]+(g-1.)*(om-rqp_io[4]))*tmp2;
		rqp_io[7]=om-cc[1];
		rqp_io[6]=rqp_io[3]*pow( (cc[1]/c[1]),gc );
		rqp_io[8]=cc[1]*cc[1]*rqp_io[6]/g;
		goto rtn_o;
	}

	/* 6 */
	if( om > ss[0][1] )	{ rqp_io[6]=rqp_io[3]; rqp_io[7]=rqp_io[4]; rqp_io[8]=rqp_io[5]; goto rtn_o; }

	rtn_o:
	
	//return alarm;
	rqp_io[6] = rqp_io[6];
	rqp_io[7] = rqp_io[7];
	rqp_io[8] = rqp_io[8];

	if(rqp_io[8] <= e_pmin)
	{
		printf("\n Final pressure <= e_pmin - no solution \n");
		printf(" Left  %e  %e  %e  \n",rqp_io[0],rqp_io[1],rqp_io[2]);
		printf(" Right %e  %e  %e  \n",rqp_io[3],rqp_io[4],rqp_io[5]);
                *alarm=79;
                return;
//		fflush(stdout);
//		exit (-1);
	}
}







